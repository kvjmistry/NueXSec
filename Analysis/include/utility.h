#ifndef UTILITY_h
#define UTILITY_h

// STD includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>

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

// Class at the top level of the selection, so most classes will be including
// from this class. Mainly provided useful functions.


class utility{

public:
    // Default Constructor
    utility();
    // ---------------------------------------------------------------------
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

    // Other definitions for code

    // For creating histogram names
    std::vector<std::string> type_prefix;

    // Cut directory names
    std::vector<std::string> cut_dirs;

    // Names of the plot types
    std::vector<std::string> plot_types;

    // Names of the classifications
    std::vector<std::string> classification_dirs;

    // enum to switch file type 
    enum type {k_mc, k_data, k_ext, k_dirt, k_variation, k_type_MAX}; 

    // enums for cut dirs
    enum enum_cut_dirs {
                k_in_fv,                           // Fiducial volume
                k_cuts_MAX
                }; 

    enum interactions {
        k_qe  = 0,
        k_res = 1,
        k_dis = 2,
        k_coh = 3,
        k_mec = 10
    };

    // enums for legend
    enum legend {
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
    
    enum pandora_classification {
        k_pandora_nu_e_other = 1,
        k_pandora_nu_e_cc0pi0p = 10,
        k_pandora_nu_e_cc0pinp = 11,
        k_pandora_nu_mu_other = 2,
        k_pandora_nu_mu_pi0 = 21,
        k_pandora_nc = 3,
        k_pandora_nc_pi0 = 31,
        k_pandora_cosmic = 4,
        k_pandora_outfv = 5,
        k_pandora_other = 6,
        k_pandora_data = 0
    };

    // enums for checking if CC or NC interaction
    enum CCNC {
        k_CC,
        k_NC
    };

}; // End Class Utility


#endif
