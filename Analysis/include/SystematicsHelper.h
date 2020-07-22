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

    // Input reweigted histogram file
    TFile *f_nuexsec;

    // Class instances
    utility _util;

    // Variables
    int run{0}, subrun{0}, event{0};
    
    TTree * tree;

    std::string run_period;

    std::vector<double> POT_v; // vector of POT for each variation 


    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char *run_period, utility _utility, const char* _mode);
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
    // Draw area norm label
    void Draw_Area_Norm(TCanvas* c);
    // -------------------------------------------------------------------------
    // Draw run period label
    void Draw_Run_Period(TCanvas* c);
    // -------------------------------------------------------------------------
    // Separate initialiser for the case of reweighted systematics.
    // We treat these different to detector systematics due to their different 
    // format compared to the detector systematics
    void InitialiseReweightingMode();
    // -------------------------------------------------------------------------
    // Draw Unisim histograms
    void PlotReweightingModeUnisim(std::string label);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------



    std::string mode{"default"}; // what mode to run this class in


    enum enum_variations {
        k_CV,
        k_bnb_diffusion,
        k_vars_MAX
    };

    enum enum_ext {
        k_NuMI,
        k_BNB,
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

    // enum for reweighter labels
    enum TH1D_xsec_label_vars {
        k_RPAup,
        k_CCMECup,
        k_AxFFCCQEup,
        k_VecFFCCQEup,
        k_DecayAngMECup,
        k_ThetaDelta2Npiup,
        k_ThetaDelta2NRadup,
        k_RPA_CCQE_Reducedup,
        k_NormCCCOHup,
        k_NormNCCOHup,
        k_RPAdn,
        k_CCMECdn,
        k_AxFFCCQEdn,
        k_VecFFCCQEdn,
        k_DecayAngMECdn,
        k_ThetaDelta2Npidn,
        k_ThetaDelta2NRaddn,
        k_RPA_CCQE_Reduceddn,
        k_NormCCCOHdn,
        k_NormNCCOHdn,
        k_weightsGenie,
        k_weightsReint,
        k_weightsPPFX,
        k_TH1D_reweighter_labels_MAX
    };

    std::vector<std::string> reweighter_labels = {
        "CV",
        "RPAup",
        "CCMECup",
        "AxFFCCQEup",
        "VecFFCCQEup",
        "DecayAngMECup",
        "ThetaDelta2Npiup",
        "ThetaDelta2NRadup",
        "RPA_CCQE_Reducedup",
        "NormCCCOHup",
        "NormNCCOHup",
        "RPAdn",
        "CCMECdn",
        "AxFFCCQEdn",
        "VecFFCCQEdn",
        "DecayAngMECdn",
        "ThetaDelta2Npidn",
        "ThetaDelta2NRaddn",
        "RPA_CCQE_Reduceddn",
        "NormCCCOHdn",
        "NormNCCOHdn",
        "weightsGenie",
        "weightsReint",
        "weightsPPFX"
    };


    // Reweighter labels for single histogram plotting
    std::vector<std::string> reweighter_labels_comb_pretty = {
        "CV",
        "RPA",
        "CC MEC",
        "Ax FF CCQE",
        "Vec FF CCQE",
        "Decay Ang MEC",
        "Theta Delta 2N pi",
        "Theta Delta 2N Rad",
        "RPA CCQE Reduced",
        "Norm CC COH",
        "Norm NC COH",
        "Genie",
        "Reinteractions",
        "PPFX"
    };
    
    std::vector<std::string> reweighter_labels_comb = {
        "CV",
        "RPA",
        "CC_MEC",
        "Ax_FF_CCQE",
        "Vec_FF_CCQE",
        "Decay_Ang_MEC",
        "Theta_Delta_2N_pi",
        "Theta_Delta_2N_Rad",
        "RPA_CCQE_Reduced",
        "Norm_CC_COH",
        "Norm_NC_COH",
        "Genie",
        "Reinteractions",
        "PPFX"
    };

    // enum for reweighter labels
    enum TH1D_xsec_label_vars_comb {
        k_comb_RPA,
        k_comb_CCMEC,
        k_comb_AxFFCCQE,
        k_comb_VecFFCCQE,
        k_comb_DecayAngMEC,
        k_comb_ThetaDelta2Npi,
        k_comb_ThetaDelta2NRad,
        k_comb_RPA_CCQE_Reduced,
        k_comb_NormCCCOH,
        k_comb_NormNCCOH,
        k_comb_weightsGenie,
        k_comb_weightsReint,
        k_comb_weightsPPFX,
        k_TH1D_reweighter_labels_comb_MAX
    };

    // enum for histogram vars
    enum TH1D_xsec_hist_vars {
        k_xsec_sel,     // Selected event histogram binned in energy
        k_xsec_bkg,     // Bkg event histogram binned in energy
        k_xsec_gen,     // Gen event histogram binned in energy
        k_xsec_sig,     // Sig event histogram binned in energy
        k_xsec_eff,     // Efficiency histogram binned in energy
        k_xsec_ext,     // EXT event histogram binned in energy
        k_xsec_dirt,    // Dirt event histogram binned in energy
        k_xsec_data,    // Data event histogram binned in energy
        k_xsec_mcxsec,  // MC Cross Section
        k_xsec_dataxsec,// Data Cross Section
        k_xsec_mcxsec_int,  // MC Cross Section Flux Integrated
        k_xsec_dataxsec_int,// Data Cross Section Flux Integrated
        k_TH1D_xsec_MAX
    };

    // Names for cross section histograms
    std::vector<std::string> xsec_types = {"sel", "bkg", "gen", "sig", "eff", "ext", "dirt", "data", "mc_xsec", "data_xsec", "mc_xsec_int", "data_xsec_int"};
    std::vector<std::string> xsec_types_pretty = {"Selected", "Background", "Generated Signal", "Signal", "Efficiency", "Beam-Off", "Dirt", "Beam-On", "MC Cross-Section", "Data Cross-Section", "MC Integrated Cross-Section", "Data Integrated Cross-Section"};

    // Containter for the central value histograms
    std::vector<TH1D*> cv_hist_vec;

    enum updn {k_up, k_dn};


}; // End Class SystematicsHelper

#endif