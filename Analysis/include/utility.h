#ifndef UTILITY_h
#define UTILITY_h

// STD includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>

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

// Class at the top level of the selection, so most classes will be including
// from this class. Mainly provided useful functions.


class utility{

public:
    // -------------------------------------------------------------------------
    // Function to configure the cut values from main.h
    std::vector<double> configure( 
                                   double _Run1_MC_POT,
                                   double _Run1_Dirt_POT,
                                   double _Run1_Data_POT,
                                   double _Run1_Data_trig,
                                   double _Run1_EXT_trig,
                                   double _Run3_MC_POT,
                                   double _Run3_Dirt_POT,
                                   double _Run3_Data_POT,
                                   double _Run3_Data_trig,
                                   double _Run3_EXT_trig,
                                   double _x1,
                                   double _x2,
                                   double _y1,
                                   double _y2,
                                   double _z1,
                                   double _z2
                                  );
    // -------------------------------------------------------------------------
    // Get a TFile from a file
    bool GetFile(TFile* &f, TString string);
    // -------------------------------------------------------------------------
    // Get a TTrees from a file
    void GetTree(TFile* f, TTree* &T, TString string);
    // -------------------------------------------------------------------------
    // Get a TDirectory from a file
    bool GetDirectory(TFile* f, TDirectory* &d, TString string);
    // -------------------------------------------------------------------------
    // Get a histogram from a file
    bool GetHist(TFile* f, TH1D* &h, TString string);
    bool GetHist(TFile* f, TH2D* &h, TString string);
    // -------------------------------------------------------------------------
    // Function to tabulate all the nuetrino types and flavours
    void Tabulate(std::string interaction, std::string classification, int type, std::vector<double> &counter_v, double weight);
    // -------------------------------------------------------------------------
    // Function to print the tabulated events
    void PrintInfo(std::vector<double> counter_v, double intime_scale_factor, double mc_scale_factor, double dirt_scale_factor, std::string cut_name, double tot_true_infv_nues, double &efficiency, double &purity);
    // -------------------------------------------------------------------------
    // Function calculate theta
    double GetTheta(double px, double py, double pz);
    // -------------------------------------------------------------------------


    // Other definitions for code
    bool verbose = false; // This should be set in the config file

    // For creating histogram names
    std::vector<std::string> type_prefix = {"MC", "Data", "EXT", "Dirt"};

    std::vector<std::string> sig_bkg_prefix = {"Signal", "Background"};

    // Cut directory names
    std::vector<std::string> cut_dirs = {
            "Unselected",    // Unselected
            "SoftwareTrig",  // Software Trigger
            "Op_Filter_PE",  // Common Optical Filter PE
            // "Op_Filter_Veto",// Common Optical Filter Michel Veto
            "Slice_ID",      // Slice ID
            "e_candidate",   // Electron Candidate
            "In_FV",         // In FV
            "Topo_Score",    // Topological Score
            "Cosmic_IP",     // Cosmic Inpact Parameter
            "Cluster_Frac",  // Cluster Fraction 
            "Shower_Score",  // Track Score
            "Shower_Contained", // Shower Containment
            "Michel_Rej",    // Michel Rejection
            "ShrHits",       // Shower Hits
            "HitRatio",      // Ratio of shr hits and slice hits
            "Moliere_Avg",   // Shower Moliere Average
            "ShrVtxDistance",// Shower to vertex distance
            "dEdx_y",          // dEdx y plane
            // "dEdx_v",          // dEdx y plane
            // "dEdx_u"           // dEdx y plane
            };


    // Names of the plot types
    std::vector<std::string> plot_types = {
                "TEff",
                "True",
                "Flash",
                "Interaction",
                "2D",
                "ParticleStack",
                "Stack"
                };

    // Names of the classifications
    std::vector<std::string> classification_dirs = {
                "nue_cc",
                "nue_cc_mixed",
                "nu_out_fv",
                "cosmic",
                "numu_cc",
                "numu_cc_pi0",
                "nc",
                "nc_pi0",
                "unmatched",
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
                k_opfilt_pe,         // Common Optical Filter PE
                // k_opfilt_veto,    // Common Optical Filter Michel Veto
                k_slice_id,          // Slice ID
                k_e_candidate,       // Electron Candidate
                k_in_fv,             // Reco Nu Vtx (SC Corr) In the FV 
                k_topo_score,        // Topo Score
                k_cosmic_ip,         // Cosmic IP
                k_cluster_frac,      // Cluster Fraction
                k_shower_score,      // Shower Score
                k_shower_contained,  // Shower containment
                k_michel_rej,        // Michel Rejection
                k_shr_hits,          // Shower Hits
                k_hit_ratio,         // Ratio of shr hits and slice hits
                k_shr_moliere_avg,   // Shower Moliere Average
                k_shr_distance,      // Shower to reco nu vertex distance
                k_dEdx_y,            // dEdx y plane
                // k_dEdx_v,            // dEdx v plane
                // k_dEdx_u,            // dEdx u plane
                k_cuts_MAX
                }; 

    // Config enums
    enum enum_config {
        k_config_Run1_MC_POT,
        k_config_Run1_Dirt_POT,
        k_config_Run1_Data_POT,
        k_config_Run1_Data_trig,
        k_config_Run1_EXT_trig,
        k_config_Run3_MC_POT,
        k_config_Run3_Dirt_POT,
        k_config_Run3_Data_POT,
        k_config_Run3_Data_trig,
        k_config_Run3_EXT_trig,
        k_config_x1,
        k_config_x2,
        k_config_y1,
        k_config_y2,
        k_config_z1,
        k_config_z2,
        k_config_MAX
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
                k_nue_cc_mixed,
                k_nu_out_fv,
                k_cosmic,
                k_numu_cc,
                k_numu_cc_pi0,
                k_nc,
                k_nc_pi0,
                k_unmatched,
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
        k_count_total_nue_cc_qe,
        k_count_total_nue_cc_res,
        k_count_total_nue_cc_dis,
        k_count_total_nue_cc_coh,
        k_count_total_nue_cc_mec,
        k_count_only_nue_cc,
        k_count_only_nue_bar_cc,
        k_count_numu_cc_qe,
        k_count_numu_cc_res,
        k_count_numu_cc_dis,
        k_count_numu_cc_coh,
        k_count_numu_cc_mec,
        k_count_tot_nue_numu_nc,
        
        k_count_nue_cc,
        k_count_nue_cc_mixed,
        k_count_nu_out_fv,
        k_count_cosmic,
        k_count_numu_cc,
        k_count_numu_cc_pi0,
        k_count_nc,
        k_count_nc_pi0,
        k_count_unmatched,
        k_count_total,
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

}; // End Class Utility


#endif
