#ifndef SYSTEMATICSHELPER_H
#define SYSTEMATICSHELPER_H

#include "Utility.h"

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
    Utility _util;

    // Variables
    int run{0}, subrun{0}, event{0};
    
    TTree * tree;

    std::vector<double> POT_v; // vector of POT for each variation 


    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise( Utility _utility);
    // -------------------------------------------------------------------------
    // Function to loop over events and calculate the cross section
    void MakeHistograms(); 
    // -------------------------------------------------------------------------
    void GetPOT();
    // -------------------------------------------------------------------------
    // Plots the variation comparisons
    void PlotVariations(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name);
    // -------------------------------------------------------------------------
    // Plots the variation comparisons
    void PlotVariationsEXT(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name);
    // -------------------------------------------------------------------------
    void SetVariationProperties(TH1D* h, int index);
    // -------------------------------------------------------------------------
    // Set options for ratio histogram
    void SetRatioOptions(TH1D* hist);
    // -------------------------------------------------------------------------
    // Draw area norm label
    void Draw_Area_Norm(TCanvas* c);
    // -------------------------------------------------------------------------
    // Separate initialiser for the case of reweighted systematics.
    // We treat these different to detector systematics due to their different 
    // format compared to the detector systematics
    void InitialiseReweightingMode();
    // -------------------------------------------------------------------------
    // Draw Unisim histograms
    void PlotReweightingModeUnisim(std::string label, int var, std::string label_pretty);
    // -------------------------------------------------------------------------
    void PlotReweightingModeMultisim(std::string label, int var, std::string label_pretty, int universes);
    // -------------------------------------------------------------------------
    // Plots the MC and data cross sections to compare them
    void CompareCVXSec(int var);
    // -------------------------------------------------------------------------
    // Function that initialises and plots the CV
    void InitialsePlotCV();
    // -------------------------------------------------------------------------
    // Compare the data cross section to mc cross section for each variation
    void CompareVariationXSec(std::string label, int var, std::string label_pretty);
    // -------------------------------------------------------------------------
    // Set the up down variation names
    void SetLabelName(std::string label, std::string &label_up, std::string &label_dn);
    // -------------------------------------------------------------------------
    // Calculate the covariance matrix for the multisims
    void CalcCovariance(std::string label, int var, std::vector<std::vector<TH1D*>> h_universe );
    // -------------------------------------------------------------------------
    // Fill the total systematic vector with the square sum of the uncertainty
    void FillSysVector(std::string variation, int var, int type, TH1D *h_up, TH1D *h_dn);
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
        k_TH1D_xsec_MAX
    };

    // enum for histogram vars
    enum TH1D_xsec_var_vars {
        k_var_integrated,     // Integrated X-Section
        k_var_reco_el_E,      // Reconstructed electron energy
        k_var_true_el_E,      // True electron energy
        // k_var_true_nu_E,      // True neutrino energy
        // k_var_reco_nu_E,      // Reconstructed neutrino energy
        k_TH1D_xsec_var_MAX
    };

    // Names for cross section histograms
    std::vector<std::string> xsec_types = {"sel", "bkg", "gen", "sig", "eff", "ext", "dirt", "data", "mc_xsec", "data_xsec"};
    std::vector<std::string> xsec_types_pretty = {"Selected", "Background", "Generated Signal", "Signal", "Efficiency", "Beam-Off", "Dirt", "Beam-On", "MC", "Data"};

    std::vector<std::string> vars = {"integrated",
                                     "reco_el_E",
                                     "true_el_E"
                                    //  "true_nu_E",
                                    //  "reco_nu_e"
                                     };


    // Use these for when we do the cross-section
    // std::vector<std::string> var_labels = {";;#nu_{e} + #bar{#nu}_{e} CC Cross-Section [10^{-39} cm^{2}]",
    //                                     ";Reco Leading Shower Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{reco}_{e}} CC Cross-Section [10^{-39} cm^{2}/GeV]",
    //                                     ";True Electron Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{true}_{e}} CC Cross-Section [10^{-39} cm^{2}/GeV]"
    //                                     // ";True #nu_{e} Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{true}_{#nu_{e}}} CC Cross-Section [10^{-39} cm^{2}/GeV]",
    //                                     // ";Reco #nu_{e} Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{reco}_{#nu_{e}}} CC Cross-Section [10^{-39} cm^{2}/GeV]"
    //                                     };
    
    // Use these for when we do the flux normalised event rate
    std::vector<std::string> var_labels = {";;#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}]",
                                        ";Reco. Leading Shower Energy [GeV];#nu_{e} + #bar{#nu}_{e} Flux Norm. Event Rate CC [cm^{2}/GeV]",
                                        ";True Electron Energy [GeV]; #nu_{e} + #bar{#nu}_{e} Flux Norm. Event Rate CC [cm^{2}/GeV]"
                                        // ";True #nu_{e} Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{true}_{#nu_{e}}} CC Cross-Section [10^{-39} cm^{2}/GeV]",
                                        // ";Reco #nu_{e} Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{reco}_{#nu_{e}}} CC Cross-Section [10^{-39} cm^{2}/GeV]"
                                        };


    std::vector<std::string> var_labels_x = {"",
                                        "Reco Leading Shower Energy [GeV]",
                                        "True Electron Energy [GeV]"
                                        // "True #nu_{e} Energy [GeV]",
                                        // "Reco #nu_{e} Energy [GeV]"
                                        };

    // Containter for the central value histograms
    std::vector<std::vector<TH1D*>> cv_hist_vec; // reco elec e, <gen, sig, etc>

    enum updn {k_up, k_dn};

    // Vectors to store the quadrature sum of the uncertainties
    // We combine all of these to then get the total error and plot it
    std::vector<std::vector<std::vector<double>>> v_genie_uni_total;   // differential variable, type, bin error [genie unisim]
    std::vector<std::vector<std::vector<double>>> v_genie_multi_total; // differential variable, type, bin error [genie multisim]
    std::vector<std::vector<std::vector<double>>> v_beamline_total;    // differential variable, type, bin error [beamline unisim]
    std::vector<std::vector<std::vector<double>>> v_hp_total;          // differential variable, type, bin error [hadron production multisim]
    std::vector<std::vector<std::vector<double>>> v_reint_total;       // differential variable, type, bin error [geant reinteraction multisim]


}; // End Class SystematicsHelper

#endif