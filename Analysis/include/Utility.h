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
    double GetCVWeight(int type, double weightSplineTimesTune, double ppfx_cv, double nu_e, int nu_pdg, bool infv, int interaction);
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
    void TurnoffIntrinsicMode(){intrinsic_mode = (char*)"empty";};
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
    // Set the axes names of the cross section plots
    void SetAxesNames(std::vector<std::string> &var_labels_xsec, std::vector<std::string> &var_labels_events,
                      std::vector<std::string> &var_labels_eff,  std::string &smear_hist_name, std::vector<std::string> &vars, double  &xsec_scale);
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
    // Apply matrix to smear from true to reco space
    void MatrixMultiply(TH1D* h_true, TH1D* &h_reco, TH2D* matrix, std::string option, bool norm);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
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
    std::vector<double> reco_shr_bins    = { 0.0, 0.30, 0.49, 0.69, 0.98, 1.47, 6.0};

    // Bins for reco shr beta
    std::vector<double> reco_shr_bins_ang = { 0.0, 6.0, 14.0, 20.5, 29.0, 42.5, 180};

    // Bins for reco shr cos beta
    std::vector<double> reco_shr_bins_cang = {-1.0, 0.6, 0.81, 0.91, 0.95, 0.98, 1.0};

    // Fine truth Binning
    std::vector<double> true_shr_bins = { 0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.15, 2.3, 2.6, 3.2, 6.0};

    std::vector<double> true_shr_bins_ang = { 0, 3, 6, 10, 12, 15, 20, 22, 24, 26, 28, 32, 36, 40, 45, 50, 60, 70, 80, 90, 100, 120, 150, 180};

    std::vector<double> true_shr_bins_cang = {-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

    // Neutrino Energy Threshold to integrate from
    // Flux bins are 0.00 ,0.06, 0.125, 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 4.00, 5.00
    double energy_threshold = 0.06; // GeV
    double elec_threshold   = 0.12; // GeV
    
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
    char * mc_file_name          = (char *)"empty";
    char * ext_file_name         = (char *)"empty";
    char * data_file_name        = (char *)"empty";
    char * dirt_file_name        = (char *)"empty";
    char * mc_file_name_out      = (char *)"empty";
    char * ext_file_name_out     = (char *)"empty";
    char * data_file_name_out    = (char *)"empty";
    char * dirt_file_name_out    = (char *)"empty";
    char * variation             = (char *)"empty";
    char * variation_file_name   = (char *)"empty";
    char * mc_tree_file_name_out = (char *)"empty";
    char * hist_file_name        = (char *)"empty";
    char * tree_file_name        = (char *)"empty";
    char * run_period            = (char *)"empty";
    char * sysmode               = (char *)"default";
    char * xsecmode              = (char *)"default";
    char * xsec_rw_mode          = (char *)"default"; // choose whether to reweight by cut or the final selection
    char * xsec_var              = (char *)"elec_E";  // What variable to do the cross section as a function of
    char * xsec_labels           = (char *)"all";
    char * uplotmode             = (char *)"default";
    char * intrinsic_mode        = (char *)"default"; // choose whether to override the nue component to accomodate the intrinsic nue sample
    char * sysplot               = (char *)"tot";     // what systematic uncertainty to plot on the CV histograms
    char * xsec_smear_mode       = (char *)"er";    // what smearing do we want to apply to the measurement? mcc8 = Marco's smearing, response = smearing using a response matrix and compare event rates
    char * xsec_bin_mode         = (char *)"standard"; // Choose whether to use standard binning (reco = truth) or fine truth binning 
    char * scale_bins            = (char *)"width"; // Choose whether to use scale the histograms by bin width. Options are standard or width
    int num_events{-1};
    int verbose{1}; // level 0 doesn't print cut summary, level 1 prints cut summary [default is 1 if unset]
    int _weight_tune{1}; // Use the GENIE Tune
    int _weight_ppfx{1}; // Use the PPFX CV Corr
    int _weight_dirt{1}; // Weight the Dirt events
    int _weight_ext{1};  // Weight the EXT events
    int _pi0_correction{1};  // The pi0 correction 0 == no correction, 1 == normalisation factor, 2 == energy dependent scaling


    // Weight configurations
    bool weight_tune{true}; // Use the GENIE Tune
    bool weight_ppfx{true}; // Use the PPFX CV Corr
    bool weight_dirt{true}; // Weight the Dirt events
    bool weight_ext{true};  // Weight the EXT events
    int  pi0_correction{1}; // The pi0 correction 0 == no correction, 1 == normalisation factor, 2 == energy dependent scaling
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
            "Unselected",     // Unselected
            "SoftwareTrig",   // Software Trigger
            "Slice_ID",       // Slice ID
            "e_candidate",    // Electron Candidate
            "In_FV",          // In FV
            "Contained_Frac", // Slice Contained Fraction
            "Topo_Score",     // Topological Score
            "Cosmic_IP",      // Pandora Cosmic Impact Parameter
            "Shower_Score",   // Track Score < 0.5
            "HitRatio",       // Ratio of shr hits and slice hits
            "Moliere_Avg",    // Shower Moliere Average
            "ShrVtxDist_dEdx_max", // 2D cut for shower to vertex distance and dedx
            "dEdx_max_no_tracks"  // dEdx all planes no tracks
            };
    
    std::vector<std::string> cut_dirs_pretty = {
            "Unselected",                                 // Unselected
            "Software Trigger",                           // Software Trigger
            "Slice ID",                                   // Slice ID
            "Electron Candidate",                         // Electron Candidate
            "In Fiducial Volume",                         // In FV
            "Contained Fraction",                         // Slice Contained Fraction
            "Topological Score",                          // Topological Score
            "Cosimc IP",                                  // Pandora Cosmic Impact Parameter
            "Shower Score",                               // Track Score < 0.5
            "Hit Ratio",                                  // Ratio of shr hits and slice hits
            "Moliere Average",                            // Shower Moliere Average
            "2D Shower Vtx Dist, dE/dx",                  // 2D cut for shower to vertex distance and dedx
            "dE/dx, 0 Tracks"                             // dEdx all planes no tracks
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
                "NC"
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
                k_swtrig,            // Software Trigger
                k_slice_id,          // Slice ID
                k_e_candidate,       // Electron Candidate
                k_in_fv,             // Reco Nu Vtx (SC Corr) In the FV 
                k_contained_frac,    // Slice Contained Fraction
                k_topo_score,        // Topo Score
                k_cosmic_ip,         // Pandora Cosmic Impact Param 3D
                k_shower_score,      // Shower Score
                k_hit_ratio,         // Ratio of shr hits and slice hits
                k_shr_moliere_avg,   // Shower Moliere Average
                k_vtx_dist_dedx,     //  2D cut for shower to vertex distance and dEdx. Only applied for > 1 track
                k_dEdx_max_no_tracks,// dEdx all planes when there is no tracks
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
        k_plot_nc,
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
        k_cut_softwaretrig,
        k_cut_nslice,
        k_cut_shower_multiplicity,
        k_cut_track_multiplicity,
        k_cut_topological_score,
        k_cut_vtx_x_sce,
        k_cut_vtx_y_sce,
        k_cut_vtx_z_sce,
        k_cut_shower_score,
        k_cut_shr_tkfit_dedx_max,
        k_cut_shr_tkfit_dedx_max_with_tracks,
        k_cut_shr_tkfit_dedx_max_no_tracks,
        k_cut_shower_to_vtx_dist,
        k_cut_hits_ratio,
        k_cut_CosmicIPAll3D,
        k_cut_contained_fraction,
        k_cut_shrmoliereavg,
        k_cut_leading_shower_theta,
        k_cut_leading_shower_phi,
        k_cut_shower_energy_cali,
        k_cut_shower_energy_cali_rebin,
        k_cut_flash_time,
        k_cut_flash_pe,
        k_cut_effective_angle,
        k_cut_effective_cosangle,
        k_cut_vars_max
    };

    // ------------------------------------------
    // variables to plot PlotVariations and SysVariations
        std::vector<std::string> vec_hist_name = {
        "h_reco_softwaretrig",
        "h_reco_nslice",
        "h_reco_shower_multiplicity",
        "h_reco_track_multiplicity",
        "h_reco_topological_score",
        "h_reco_vtx_x_sce",
        "h_reco_vtx_y_sce",
        "h_reco_vtx_z_sce",
        "h_reco_shower_score",
        "h_reco_shr_tkfit_dedx_max",
        "h_reco_shr_tkfit_dedx_max_with_tracks",
        "h_reco_shr_tkfit_dedx_max_no_tracks",
        "h_reco_shower_to_vtx_dist",
        "h_reco_hits_ratio",
        "h_reco_CosmicIPAll3D",
        "h_reco_contained_fraction",
        "h_reco_shrmoliereavg",
        "h_reco_leading_shower_theta",
        "h_reco_leading_shower_phi",
        "h_reco_shower_energy_cali",
        "h_reco_shower_energy_cali_rebin",
        "h_reco_flash_time",
        "h_reco_flash_pe",
        "h_reco_effective_angle",
        "h_reco_effective_cosangle"
    };
 
    // x axis label for those plots
    std::vector<std::string> vec_axis_label = {
        "Software Trigger",
        "Pandora Slice ID",
        "Shower Multiplicity",
        "Track Multiplicity",
        "Topological Score",
        "Reco Vertex X [cm]",
        "Reco Vertex Y [cm]",
        "Reco Vertex Z [cm]",
        "Shower Score",
        "Leading Shower dE/dx (All Planes) [MeV/cm]",
        "Leading Shower dE/dx (All Planes) (with tracks) [MeV/cm]",
        "Leading Shower dE/dx (All Planes) (0 tracks) [MeV/cm]",
        "Leading Shower to Vertex Distance [cm]",
        "Hit Ratio",
        "Pandora Cosmic Impact Parameter 3D [cm]",
        "Contained Fraction (PFP hits in FV / hits in slice)",
        "Leading Shower Moliere Average [deg]",
        "Leading Shower Theta [deg]",
        "Leading Shower Phi [deg]",
        "Reconstructed Leading Shower Energy [GeV]",
        "Reconstructed Leading Shower Energy [GeV]",
        "Largest Flash Time [#mus]",
        "Largest Flash Intensity [PE]",
        "Leading Shower Effective Angle [deg]",
        "Leading Shower Cosine Effective Angle [deg]"
    };

    // list of detector variations
    std::vector<std::string> vec_var_string = {
        "CV",
        //"BNB_Diffusion",
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

}; // End Class Utility


#endif
