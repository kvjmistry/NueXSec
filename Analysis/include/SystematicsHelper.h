#ifndef SYSTEMATICSHELPER_H
#define SYSTEMATICSHELPER_H

#include "Utility.h"
#include "WienerSVD.h"

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

    WienerSVD _wSVD;

    std::vector<double> POT_v; // vector of POT for each variation 

    bool scale_bins = false;

    double Data_POT;
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
    void PlotVariationsEXT(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name);
    // -------------------------------------------------------------------------
    // Plots the Sys Variations
    void SysVariations(int hist_index, const char* print_name, int cut, const char* x_axis_name);
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
    // This function will get the reweighted plots by cut and save the uncertainties to a file
    // so we can plot the total systematic uncertainty
    void InitialiseReweightingModeCut();
    // -------------------------------------------------------------------------
    // Draw Unisim histograms
    void PlotReweightingModeUnisim(std::string label, int var, std::string label_pretty);
    // -------------------------------------------------------------------------
    void PlotReweightingModeMultisim(std::string label, int var, std::string label_pretty, int universes);
    // -------------------------------------------------------------------------
    // Plots the MC and data cross sections to compare them
    void CompareCVXSec();
    // -------------------------------------------------------------------------
    // Plots the MC and data cross sections to compare them
    // This is with no ratio plot
    void CompareCVXSecNoRatio();
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
    // Calculate the covariance, correlation and fractional covariance matrices
    void CalcMatrices(std::string label, int var, std::vector<std::vector<TH1D*>> h_universe, int _type, TH1D* h_CV  );
    // -------------------------------------------------------------------------
    // Fill the total systematic vector with the square sum of the uncertainty
    void FillSysVector(std::string variation, int var, int type, TH1D *h_up, TH1D *h_dn, TH1D* h_CV);
    // -------------------------------------------------------------------------
    // Fill vector with the statistical uncertainties
    void FillStatVector();
    // -------------------------------------------------------------------------
    // Fill POT counting uncertainty vector
    void FillPOTCountingVector();
    // -------------------------------------------------------------------------
    // Print Summary of uncertainties
    void PrintUncertaintySummary();
    // -------------------------------------------------------------------------
    // Resize the containers for storing the final uncertainties
    void InitialiseUncertaintyVectors();
    // -------------------------------------------------------------------------
    // Save the covariance matrix
    void SaveCovMatrix(TH2D* cov, std::string print_name);
    // -------------------------------------------------------------------------
    // Make the total beamline sys error plots
    void PlotTotUnisim(std::string unisim_type);
    // -------------------------------------------------------------------------
    // Set the fill colours of the unisim variations
    void SetUnisimColours(std::string label, TH1D* h_up, TH1D* h_dn);
    // -------------------------------------------------------------------------
    // Get the systematic uncertainty for each cut for a specific set of variables
    void GetCutSysUncertainty(std::string histname, int cut_index, std::string label, int num_uni, std::string var_type, TH1D* &h_err);
    // -------------------------------------------------------------------------
    // Plot detector variation histograms for the cross section variables
    void PlotReweightingModeDetVar(std::string label, int var, int detvar_index, std::string label_pretty);
    // -------------------------------------------------------------------------
    // Write the cut histograms to file
    void SaveCutHistograms(std::vector<std::tuple<std::string, int, std::string>> tuple_label);
    // -------------------------------------------------------------------------
    // Write the cut histograms to file for detector variations
    void SaveCutHistogramsDetVar();
    // -------------------------------------------------------------------------
    // Make a plot of the systematics in one plot
    void MakeTotUncertaintyPlot();
    // -------------------------------------------------------------------------
    // Initialse the matrix of covariance matrices
    void InitialseCovarianceVector();
    // -------------------------------------------------------------------------
    // Write the results to a file for ease of storage
    void ExportResult(TFile* f);
    // -------------------------------------------------------------------------
    // Save the total data cross section result to file
    void ExportTotalCrossSectionResult();
    // -------------------------------------------------------------------------
    // Change units of covariane matrix so we can add them
    void ConvertCovarianceUnits(TH2D* &h_cov, TH1D *h_input, TH1D* h_output);
    // -------------------------------------------------------------------------



    std::string mode{"default"}; // what mode to run this class in


    enum enum_variations {
        k_CV,
        k_LYRayleigh,
        k_LYDown,
        // k_LYAttenuation,
        k_SCE,
        k_Recomb2,
        k_WireModX,
        k_WireModYZ,
        k_WireModThetaXZ,
        k_WireModThetaYZ_withSigmaSplines,
        // k_WireModThetaYZ_withoutSigmaSplines,
        // k_WireModdEdX,
        k_vars_MAX
    };

    enum enum_ext {
        k_NuMI,
        k_BNB,
        k_ext_MAX
    };

    std::vector<std::string> var_string = {
        "CV",
        "LYRayleigh",
        "LYDown",
        // "LYAttenuation",
        "SCE",
        "Recomb2",
        "WireModX",
        "WireModYZ",
        "WireModThetaXZ",
        "WireModThetaYZ_withSigmaSplines",
        // "WireModThetaYZ_withoutSigmaSplines",
        // "WireModdEdX"
    };

    // strings used for the legend of the plots of the detector variations
    std::vector<std::string> var_string_pretty = {
        "CV",
        "LY Rayleigh",
        "LY Down",
        // "LY Attenuation",
        "SCE",
        "Recombination",
        "WM X",
        "WM YZ",
        "WM Theta XZ",
        "WM Theta YZ w/ Spl.",
        // "WM Theta YZ w/o Spl.",
        // "WM dE/dX" 
    };

    // colors used for the plots of the detector variations
    std::vector<int> var_string_pretty_color = {
         1,  // CV
         2,  // LYRayleigh
         95, // LY Down
         6,  // LYAttenuation
         3,  // SCE
         9,  // Recombination
         30, // WM X
         4,  // WM YZ
         53, // WM Theta XZ
         8,  // WM Theta YZ w/ Spl.
        //  8, // WM Theta YZ w/o Spl.
        //  7 // WM dE/dX
    };

    // enum for histogram types
    enum TH1D_xsec_hist_vars {
        k_xsec_sel,          // Selected event histogram binned in energy
        k_xsec_bkg,          // Bkg event histogram binned in energy
        k_xsec_gen,          // Gen event histogram binned in energy
        k_xsec_gen_smear,    // Gen event histogram binned in energy with smeared truth
        k_xsec_sig,          // Sig event histogram binned in energy
        k_xsec_eff,          // Efficiency histogram binned in energy
        k_xsec_ext,          // EXT event histogram binned in energy
        k_xsec_dirt,         // Dirt event histogram binned in energy
        k_xsec_data,         // Data event histogram binned in energy
        k_xsec_mcxsec,       // MC Cross Section
        k_xsec_mcxsec_smear, // MC Cross Section smeared truth
        k_xsec_dataxsec,     // Data Cross Section
        k_TH1D_xsec_MAX
    };

    // enum for histogram vars
    enum TH1D_xsec_var_vars {
        k_var_integrated,     // Integrated X-Section
        k_var_reco_el_E,      // Reconstructed electron energy
        k_var_true_el_E,      // True electron energy
        k_TH1D_xsec_var_MAX
    };

    // Names for cross section histograms
    std::vector<std::string> xsec_types = {"sel", "bkg", "gen", "gen_smear", "sig", "eff", "ext", "dirt", "data", "mc_xsec", "mc_xsec_smear", "data_xsec"};
    std::vector<std::string> xsec_types_pretty = {"Selected", "Background", "Generated Signal", "Smeared Prediction", "Signal", "Efficiency", "Beam-Off", "Dirt", "Beam-On", "MC", "MC Smear",  "Data"};

    std::vector<std::string> vars = {"integrated","recoX", "trueX" };

    // Choose the cross section scale to set the histogram
    double xsec_scale = 13.0;
    
    // Use these for when we do the flux normalised event rate
    // std::vector<std::string> var_labels_xsec = {";;#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}]",
    //                                        ";Reco. Leading Shower Energy [GeV];#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}/GeV]",
    //                                        ";True e#lower[-0.5]{-} + e^{+} Energy [GeV];#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}/GeV]"
                                        // };
    
    std::vector<std::string> var_labels_xsec = {};

    std::vector<std::string> var_labels_events = {};

    std::vector<std::string> var_labels_eff = {};

    std::string smear_hist_name = ";True e#lower[-0.5]{-} + e^{+} Energy [GeV];Leading Shower Energy [GeV]";

    // Containter for the central value histograms
    std::vector<std::vector<TH1D*>> cv_hist_vec; // reco elec e, <gen, sig, etc>

    enum updn {k_up, k_dn};

    // Names of the overall errors
    enum enum_sys {
        k_err_tot,
        k_err_stat,
        k_err_sys,
        k_err_genie_uni,
        k_err_genie_multi,
        k_err_hp,
        k_err_beamline,
        k_err_dirt,
        k_err_pot,
        k_err_reint,
        k_err_detvar,
        k_err_pi0,
        k_err_mcstats,
        k_ERR_MAX
    };

    std::vector<std::string> systematic_names = {
        "tot",
        "stat",
        "sys",
        "genie_uni",
        "genie_multi",
        "hp",
        "beamline",
        "dirt",
        "pot",
        "reint",
        "detvar",
        "pi0",
        "mcstats"
    };

    // Vector to store the quadrature sum of the uncertainties
    std::vector<std::vector<std::vector<std::vector<double>>>> v_err; // error type -- var -- type -- bin error . Units are percent squared. To ger the raw err, sqrt and * 0.01 then multiply by the bin content
    
    // Vector to store all the final covariance matrices
    std::vector<std::vector<std::vector<TH2D*>>> h_cov_v; // var -- type -- label


    // Vector to store total uncertainty histograms
    std::vector<std::vector<std::vector<TH1D*>>> h_cut_err; // label -- cut -- variable


}; // End Class SystematicsHelper

#endif
