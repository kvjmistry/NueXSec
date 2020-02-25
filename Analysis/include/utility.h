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

// Class at the top level of the selection, so most classes will be including
// from this class. Mainly provided useful functions.


class utility{

public:
    // -------------------------------------------------------------------------
    // Function to configure the cut values from main.h
    std::vector<double> configure_cuts( double _x1,
                                        double _x2,
                                        double _y1,
                                        double _y2,
                                        double _z1,
                                        double _z2,
                                        double flash_pe_threshold,
                                        double flash_time_start,
                                        double flash_time_end,
                                        double tolerance,
                                        double shwr_nue_tolerance,
                                        double trk_nue_tolerance,
                                        double shwr_hit_threshold,
                                        double shwr_hit_threshold_collection,
                                        double tolerance_open_angle_min,
                                        double tolerance_open_angle_max,
                                        double tolerance_dedx_min,
                                        double tolerance_dedx_max,
                                        double dist_tolerance,
                                        double pfp_hits_length_tolerance,
                                        double ratio_tolerance,
                                        bool do_variations
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
    // -------------------------------------------------------------------------
    // Function to tabulate all the nuetrino types and flavours
    void Tabulate(std::string interaction, std::string classification, int type, std::vector<int> &counter_v);
    // -------------------------------------------------------------------------
    // Function to print the tabulated events
    void PrintInfo(std::vector<int> counter_v, double intime_scale_factor, double data_scale_factor, double dirt_scale_factor, std::string cut_name, int tot_true_infv_nues);
    // -------------------------------------------------------------------------

    // Other definitions for code
    bool verbose = false; // This should be set in the config file

    // For creating histogram names
    std::vector<std::string> type_prefix = {"MC", "Data", "EXT", "Dirt"};

    // Cut directory names
    std::vector<std::string> cut_dirs = {
            "Unselected",  // Unselected
            "Slice_ID",    // Slice ID
            };


    // Names of the plot types
    std::vector<std::string> plot_types = {
                "TEff",
                "2D",
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

    // enum to switch file type 
    enum type {k_mc, k_data, k_ext, k_dirt, k_variation, k_type_MAX}; 

    // enums for cut dirs
    enum enum_cut_dirs {
                k_unselected, // Unselected 
                k_slice_id,      // Slice ID
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

}; // End Class Utility


#endif
