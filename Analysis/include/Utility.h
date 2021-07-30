#ifndef UTILITY_H
#define UTILITY_H

// STD includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <iomanip>

// Root Includes
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TCut.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TColor.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TChain.h"
#include "TDecompSVD.h"
#include "TObjString.h"

/*

Class at the top level of the selection, so most classes will be including
from this class. Mainly provided useful functions which are common to all classes.
Main idea is that we dont repeat functions in multiple places

*/


class Utility{

public:
    // -------------------------------------------------------------------------
    // Initalise variables
    void Initalise(int argc, char *argv[], std::string usage,std::string usage2, std::string usage3);
    // -------------------------------------------------------------------------
    // Get a TFile from a file
    bool GetFile(TFile* &f, TString string);
    // -------------------------------------------------------------------------
    // Get a TTrees from a file
    bool GetTree(TFile* f, TTree* &T, TString string);
    // -------------------------------------------------------------------------
    // Get a TDirectory from a file
    bool GetDirectory(TFile* f, TDirectory* &d, TString string);
    // -------------------------------------------------------------------------
    // Get a 1D histogram from a file
    bool GetHist(TFile* f, TH1D* &h, TString string);
    // -------------------------------------------------------------------------
    // Get a 2D histogram from a file
    bool GetHist(TFile* f, TH2D* &h, TString string);
    // -------------------------------------------------------------------------
    // Check whether a weight has a suitable value, templated with double and float
    template<typename T> void CheckWeight(T &weight);
    // -------------------------------------------------------------------------
    // Get the CV weight correction
    double GetCVWeight(int type, double weightSplineTimesTune, double ppfx_cv, double nu_e, int nu_pdg, bool infv, int interaction, double elec_e);
    // -------------------------------------------------------------------------
    // Get the pi0 weight correction
    void GetPiZeroWeight(double &weight, int pizero_mode, int nu_pdg, int ccnc, int npi0, double pi0_e);
    // -------------------------------------------------------------------------
    // Create another directory in the plots folder
    void CreateDirectory(std::string folder);
    // -------------------------------------------------------------------------
    // Function to tabulate all the nuetrino types and flavours
    void Tabulate(bool inFV, std::string interaction, std::string classification, std::string pi0_classification, int type, std::vector<double> &counter_v, double weight);
    // -------------------------------------------------------------------------
    // Function calculate theta
    double GetNuMIAngle(double px, double py, double pz, std::string direction);
    // -------------------------------------------------------------------------
    // Check if vertex is in the FV
    bool in_fv(double x, double y, double z);
    // -------------------------------------------------------------------------
    // Increase the label size of 1D histograms
    void IncreaseLabelSize(TH1D *h, TCanvas *c);
    // -------------------------------------------------------------------------
    // Increase the label size of 2D histograms
    void IncreaseLabelSize(TH2D *h, TCanvas *c);
    // -------------------------------------------------------------------------
    // Draw the run period on the plot
    void Draw_Run_Period(TCanvas *c, double x1, double y1, double x2, double y2);
    // -------------------------------------------------------------------------
    // Draw the data to MC ratio on the plot
    void Draw_Data_MC_Ratio(TCanvas *c, double ratio, double x1, double y1, double x2, double y2);
    // -------------------------------------------------------------------------
    // Draw the Data POT on the plot
    void Draw_Data_POT(TCanvas *c, double pot, double x1, double y1, double x2, double y2);
    // -------------------------------------------------------------------------
    // Draw MicroBooNE Simulation on canvas
    void Draw_ubooneSim(TCanvas *c, double x1, double y1, double x2, double y2);
    // -------------------------------------------------------------------------
    // Draw FHC/RHC on the plots
    void Draw_Nu_Mode(TCanvas* c, double x1, double y1, double x2, double y2);
    // -------------------------------------------------------------------------
    // Function to customise the TLatex
    void SetTextProperties(TLatex* text);
    // -------------------------------------------------------------------------
    // Initialise the TPad size for a ratio type of plot
    void SetTPadOptions(TPad *topPad, TPad *bottomPad);
    // -------------------------------------------------------------------------
    // Check if POT from input file matches with the value in config.txt
    void CheckPOT();
    // -------------------------------------------------------------------------
    // Check if a specific histogram exists in the given vector of strings
    bool CheckHistogram(std::vector<std::string> vector, TString hist_name);
    // -------------------------------------------------------------------------
    // Turn off the intrinsic nue mode
    void TurnoffIntrinsicMode(){intrinsic_mode = (char*)"default";};
    // -------------------------------------------------------------------------
    // Calculate a covariance matrix
    void CalcCovariance(std::vector<TH1D*> h_universe, TH1D *h_CV, TH2D *h_cov);
    // -------------------------------------------------------------------------
    // Calculate a correlation matrix
    void CalcCorrelation(TH1D *h_CV, TH2D  *h_cov, TH2D *h_cor);
    // -------------------------------------------------------------------------
    // Calculate a fractional covariance matrix
    // ** Requires the input frac cov to be a cloned version of the covariance matrix **
    void CalcCFracCovariance(TH1D *h_CV, TH2D *h_frac_cov);
    // -------------------------------------------------------------------------
    // Calulate a chi squared using a covarinace matrix, for a model to data
    void CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval);
    // -------------------------------------------------------------------------
    // Calculate a chi squared using a covarinace matrix, for a model to data with no correlations between bins
    void CalcChiSquaredNoCorr(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval);
    // -------------------------------------------------------------------------
    // Calculate the chi-quared with a list of bin indexes to remove
    void CalcChiSquaredRemove(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval, std::vector<int> indexes);
    // -------------------------------------------------------------------------
    // Set the axes names of the cross section plots
    void SetAxesNames();
    // -------------------------------------------------------------------------
    // Undo a bin width scaling
    void UndoBinWidthScaling(TH1D* &hist);
    // -------------------------------------------------------------------------
    // Save a 1D hsitogram as a PDF
    void Save1DHists(const char *print_name, TH1D* hist, const char* draw_option);
    // -------------------------------------------------------------------------
    // Save a 2D hsitogram as a PDF
    void Save2DHists(const char *print_name, TH2D* hist, const char* draw_option);
    // -------------------------------------------------------------------------
    // Convert 2D histogram to bin index and then save it to pdf
    void Save2DHistsBinIndex(const char *print_name, TH2D* hist, const char* draw_option, std::string type);
    // -------------------------------------------------------------------------
    // Apply matrix to smear from true to reco space
    void MatrixMultiply(TH1D* h_true, TH1D* &h_reco, TH2D* matrix, std::string option, bool norm);
    // -------------------------------------------------------------------------
    // Change units of covariane matrix so we can add them
    void ConvertCovarianceUnits(TH2D* &h_cov, TH1D *h_input, TH1D* h_output);
    // -------------------------------------------------------------------------
    // Convert a covariance matrix from un-bin width normalised to bin-width normalised units
    void ConvertCovarianceBinWidth(TH2D* &h_cov, TH1D *h_input);
    // -------------------------------------------------------------------------


    // Variables
    std::string red     = "\033[0;31m";
    std::string green   = "\033[0;32m";
    std::string yellow  = "\033[0;33m";
    std::string blue    = "\033[0;34m";
    std::string magenta = "\033[0;35m";
    std::string cyan    = "\033[0;36m";
    std::string reset   = "\033[0m";

    // Bins for the reconstructed shower energy
    std::vector<double> reco_shr_bins    = { 0.0, 0.30, 0.47, 0.70, 0.99, 1.43, 3.0, 6.0};

    // Bins for reco shr beta
    std::vector<double> reco_shr_bins_ang = { 0.0, 6.0, 13.5, 20.0, 27.5, 39.5, 180};

    // Bins for reco shr cos beta
    std::vector<double> reco_shr_bins_cang = {-1.0, 0.6, 0.79, 0.90, 0.95, 1.0};
    // std::vector<double> reco_shr_bins_cang = {-1.0, 0.6, 0.81, 0.91, 0.95, 1.0};

    // Fine truth Binning
    std::vector<double> true_shr_bins = { 0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.15, 2.3, 2.6, 3.2, 6.0};

    std::vector<double> true_shr_bins_ang = { 0, 3, 6, 10, 12, 15, 20, 22, 24, 26, 28, 32, 36, 40, 45, 50, 60, 70, 80, 90, 100, 120, 150, 180};

    std::vector<double> true_shr_bins_cang = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

    // Neutrino Energy Threshold to integrate from
    // Flux bins are 0.00 ,0.06, 0.125, 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 4.00, 5.00
    double energy_threshold = 0.06; // GeV
    double elec_threshold   = 0.12; // GeV

    // For MCC8 values
    // double energy_threshold = 0.25; // GeV
    // double elec_threshold   = 0.0; // GeV
    
    
    bool slim                      = false;
    bool make_histos               = false;
    bool run_selection             = true;
    bool area_norm                 = false;
    bool calc_cross_sec            = false;
    bool overwritePOT              = false; 
    bool run_sys                   = false;
    bool run_uplot                 = false;
    bool print                     = false;
    bool print_mc                  = false;
    bool print_data                = false;
    bool print_ext                 = false;
    bool print_dirt                = false;
    bool plot_sys_uncertainty      = false;
    bool use_gpvm                  = false; // choose whether you are on the gpvm or Krish's personal laptop

    // inputs 
    char * mc_file_name            = (char *)"empty";
    char * ext_file_name           = (char *)"empty";
    char * data_file_name          = (char *)"empty";
    char * dirt_file_name          = (char *)"empty";
    char * fake_intrinsic_file     = (char *)"empty";
    char * mc_file_name_out        = (char *)"empty";
    char * ext_file_name_out       = (char *)"empty";
    char * data_file_name_out      = (char *)"empty";
    char * dirt_file_name_out      = (char *)"empty";
    char * variation               = (char *)"empty";
    char * variation_file_name     = (char *)"empty";
    char * mc_tree_file_name_out   = (char *)"empty";
    char * ext_tree_file_name_out  = (char *)"empty";
    char * data_tree_file_name_out = (char *)"empty";
    char * dirt_tree_file_name_out = (char *)"empty";
    char * hist_file_name          = (char *)"empty";
    char * tree_file_name          = (char *)"empty";
    char * run_period              = (char *)"empty";
    char * sysmode                 = (char *)"default";
    char * xsecmode                = (char *)"default";
    char * xsec_rw_mode            = (char *)"default";  // choose whether to reweight by cut or the final selection
    char * xsec_var                = (char *)"elec_E";   // What variable to do the cross section as a function of
    char * xsec_labels             = (char *)"all";
    char * uplotmode               = (char *)"default";
    char * intrinsic_mode          = (char *)"default";  // choose whether to override the nue component to accomodate the intrinsic nue sample
    char * sysplot                 = (char *)"tot";      // what systematic uncertainty to plot on the CV histograms
    char * xsec_smear_mode         = (char *)"er";       // what smearing do we want to apply to the measurement? mcc8 = Marco's smearing, response = smearing using a response matrix and compare event rates
    char * xsec_bin_mode           = (char *)"standard"; // Choose whether to use standard binning (reco = truth) or fine truth binning 
    char * scale_bins              = (char *)"standard";    // Choose whether to use scale the histograms by bin width. Options are standard or width
    char * fakedataname            = (char *)"empty";
    int num_events{-1};
    int verbose{1}; // level 0 doesn't print cut summary, level 1 prints cut summary [default is 1 if unset]
    int _weight_tune{1}; // Use the GENIE Tune
    int _weight_ppfx{0}; // Use the PPFX CV Corr
    int _weight_dirt{1}; // Weight the Dirt events
    int _weight_ext{1};  // Weight the EXT events
    int _pi0_correction{2};  // The pi0 correction 0 == no correction, 1 == normalisation factor, 2 == energy dependent scaling

    bool zoom{false};        // bool to decide whether to zoom in on the plots
    bool isfakedata{false};  // bool for using MC as fake data
    bool isvariation{false}; // Is some different input to the selection
    bool usefluggflux{false};


    // Weight configurations
    bool weight_tune{true}; // Use the GENIE Tune
    bool weight_ppfx{false}; // Use the PPFX CV Corr
    bool weight_dirt{true}; // Weight the Dirt events
    bool weight_ext{true};  // Weight the EXT events
    int  pi0_correction{2}; // The pi0 correction 0 == no correction, 1 == normalisation factor, 2 == energy dependent scaling
    bool tune_mec{false};

    // Scale factors to scale samples to data)
    double mc_scale_factor     = 1.0;
    double ext_scale_factor    = 1.0;
    double dirt_scale_factor   = 1.0;
    double intrinsic_weight    = 1.0;

    // POT 
    std::vector<double> config_v;

    std::vector<std::string> confignames = {
        "Run1_MC_POT",
        "Run1_Dirt_POT",
        "Run1_Data_POT",
        "Run1_Data_trig",
        "Run1_EXT_trig",
        "Run1_Intrinsic_POT",
        "Run3_MC_POT",
        "Run3_Dirt_POT",
        "Run3_Data_POT",
        "Run3_Data_trig",
        "Run3_EXT_trig",
        "Run3_Intrinsic_POT",
        "x1", "x2", "y1", "y2", "z1", "z2"
    };

    // enums for config variables
    enum enum_config {
                k_Run1_MC_POT,
                k_Run1_Dirt_POT,
                k_Run1_Data_POT,
                k_Run1_Data_trig,
                k_Run1_EXT_trig,
                k_Run1_Intrinsic_POT,
                k_Run3_MC_POT,
                k_Run3_Dirt_POT,
                k_Run3_Data_POT,
                k_Run3_Data_trig,
                k_Run3_EXT_trig,
                k_Run3_Intrinsic_POT,
                k_config_x1,
                k_config_x2,
                k_config_y1,
                k_config_y2,
                k_config_z1,
                k_config_z2,
                k_config_MAX
                };

    // For creating histogram names
    std::vector<std::string> type_prefix = {"MC", "Data", "EXT", "Dirt"};

    std::vector<std::string> sig_bkg_prefix = {"Signal", "Background"};

    // Cut directory names
    std::vector<std::string> cut_dirs = {
            "Unselected"     // Unselected
            };
    
    std::vector<std::string> cut_dirs_pretty = {
            "Unselected"                                 // Unselected
            };


    // Names of the plot types
    std::vector<std::string> plot_types = {
                "TEff",
                "True",
                "Flash",
                "Interaction",
                "2D",
                "pizero",
                "numu",
                "ParticleStack",
                "Stack"
                };

    // Names of the classifications
    std::vector<std::string> classification_dirs = {
                "nue_cc",
                "nuebar_cc",
                "nu_out_fv",
                "cosmic",
                "cosmic_nue", 
                "cosmic_nuebar",
                "numu_cc",
                "numu_cc_pi0",
                "nc",
                "nc_pi0",
                "unmatched",
                "unmatched_nue",
                "unmatched_nuebar",
                "thr_nue",    // below threshold nue 
                "thr_nuebar", // below threshold nuebar  
                "ext",
                "data",
                "dirt"
                };

    // Names of the genie interaction types
    std::vector<std::string> interaction_types = {
                "CC_QE",
                "CC_Res",
                "CC_DIS",
                "CC_Coh",
                "CC_MEC",
                "CC_Tot"
                };
    
     // Names of the Particle types
    std::vector<std::string> particle_types = {
                "electron",
                "muon",
                "proton",
                "photon",
                "kaon",
                "pion",
                "cosmic",
                "neutron",
                "unmatched",
                "ext",
                "data",
                "dirt"
                };

    // enum to switch file type 
    enum type {k_mc, k_data, k_ext, k_dirt, k_type_MAX}; 

    // enums for cut dirs
    enum enum_cut_dirs {
                k_unselected,        // Unselected 
                k_cuts_MAX
                }; 

    // Genie interaction enums
    enum enum_interactions {
        k_qe  = 0,
        k_res = 1,
        k_dis = 2,
        k_coh = 3,
        k_mec = 10
    };

    // Genie interaction enums
    enum enum_plot_interactions {
        k_plot_qe,
        k_plot_res,
        k_plot_dis,
        k_plot_coh,
        k_plot_mec,
        k_plot_tot,
        k_interactions_MAX
    };

    // enums for legend
    enum enum_classification {
                k_nue_cc,
                k_nuebar_cc,
                k_nu_out_fv,
                k_cosmic,
                k_cosmic_nue,
                k_cosmic_nuebar,
                k_numu_cc,
                k_numu_cc_pi0,
                k_nc,
                k_nc_pi0,
                k_unmatched,
                k_unmatched_nue,
                k_unmatched_nuebar,
                k_thr_nue,
                k_thr_nuebar,
                k_leg_ext,
                k_leg_data,
                k_leg_dirt,
                k_classifications_MAX
                };
    
    // The PeLEE teams classifciation of events
    enum enum_pandora_classification {
        k_pandora_nu_e_other   = 1,
        k_pandora_nu_e_cc0pi0p = 10,
        k_pandora_nu_e_cc0pinp = 11,
        k_pandora_nu_mu_other  = 2,
        k_pandora_nu_mu_pi0    = 21,
        k_pandora_nc           = 3,
        k_pandora_nc_pi0       = 31,
        k_pandora_cosmic       = 4,
        k_pandora_outfv        = 5,
        k_pandora_other        = 6,
        k_pandora_data         = 0
    };

    // enums for checking if CC or NC interaction
    enum enum_CCNC {
        k_CC,
        k_NC
    };

    // enum for counter vector
    enum enum_counters {
        k_count_nue_cc_qe,
        k_count_nue_cc_res,
        k_count_nue_cc_dis,
        k_count_nue_cc_coh,
        k_count_nue_cc_mec,
        
        k_count_nuebar_cc_qe,
        k_count_nuebar_cc_res,
        k_count_nuebar_cc_dis,
        k_count_nuebar_cc_coh,
        k_count_nuebar_cc_mec,
        
        k_count_nue_cc_infv,
        k_count_nuebar_cc_infv,
        k_count_nue_cc_incryo,
        k_count_nuebar_cc_incryo,
        
        k_count_numu_cc_qe,
        k_count_numu_cc_res,
        k_count_numu_cc_dis,
        k_count_numu_cc_coh,
        k_count_numu_cc_mec,
        
        k_count_numubar_cc_qe,
        k_count_numubar_cc_res,
        k_count_numubar_cc_dis,
        k_count_numubar_cc_coh,
        k_count_numubar_cc_mec,
        
        k_count_numu_cc_infv,
        k_count_numubar_cc_infv,
        k_count_numu_cc_incryo,
        k_count_numubar_cc_incryo,

        k_count_pi0_nue_cc_nopi0,
        k_count_pi0_nue_cc_pi0,
        k_count_pi0_nuebar_cc_nopi0,
        k_count_pi0_nuebar_cc_pi0,
        k_count_pi0_numu_cc_nopi0,
        k_count_pi0_numu_cc_pi0,
        k_count_pi0_nc_nopi0,
        k_count_pi0_nc_pi0,
        
        k_count_nue_cc,
        k_count_nuebar_cc,
        k_count_nu_out_fv,
        k_count_cosmic,
        k_count_numu_cc,
        k_count_numu_cc_pi0,
        k_count_nc,
        k_count_nc_pi0,
        k_count_unmatched,
        k_count_unmatched_nue,
        k_count_cosmic_nue,
        k_count_unmatched_nuebar,
        k_count_cosmic_nuebar,
        k_count_total_mc,
        k_count_thr_nue,
        k_count_thr_nuebar,
        
        k_count_data,
        k_count_ext,
        k_count_dirt,
        k_COUNTER_MAX
    };

    // Enums for singal and background separation plots
    enum enum_sig_bkg {
        k_signal,
        k_background,
        k_sig_bkg_MAX
    };

    // Enum for labelling the particles by type
    enum enum_particle_type {
        k_electron,
        k_muon,
        k_proton,
        k_photon,
        k_kaon,
        k_pion,
        k_part_cosmic,
        k_neutron,
        k_part_unmatched,
        k_part_ext,
        k_part_data,
        k_part_dirt,
        k_particles_MAX
    };

    // Variables to get the systematic uncertainty by cut for
    enum TH1D_cut_vars {
        k_pi0_mass,
        k_cut_vars_max
    };

    // ------------------------------------------
    // variables to plot PlotVariations and SysVariations
        std::vector<std::string> vec_hist_name = {
        "h_pi0_mass"
    };
 
    // x axis label for those plots
    std::vector<std::string> vec_axis_label = {
        "Pi0 Mass [MeV]"
    };

    // list of detector variations
    std::vector<std::string> vec_var_string = {
        "CV",
        "LYRayleigh",
        "LYAttenuation",
        "SCE",
        "Recomb2",
        "WireModX",
        "WireModYZ",
        "WireModThetaXZ",
        "WireModThetaYZ_withSigmaSplines",
        "WireModThetaYZ_withoutSigmaSplines",
        "WireModdEdX"
   };


    // Global access to histogram names
    std::vector<std::string> var_labels_xsec = {};

    std::vector<std::string> var_labels_events = {};

    std::vector<std::string> var_labels_eff = {};

    std::string smear_hist_name;
    
    std::string ac_hist_name;

    std::vector<std::string> vars = {};

    double xsec_scale;

}; // End Class Utility


#endif
