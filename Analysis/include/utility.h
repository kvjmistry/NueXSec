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

    // ---------------------------------------------------------------------
    // Get a TFile from a file
    bool GetFile(TFile* &f, TString string);
    // ---------------------------------------------------------------------
    // Get a TTrees from a file
    void GetTree(TFile* f, TTree* &T, TString string);
    // ---------------------------------------------------------------------
    // Get a TDirectory from a file
    bool GetDirectory(TFile* f, TDirectory* &d, TString string);
    // ---------------------------------------------------------------------
    // Get a histogram from a file
    bool GetHist(TFile* f, TH1D* &h, TString string);
    // ---------------------------------------------------------------------

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

    // enums for legend
    enum legend {
                k_nue_cc,
                k_nue_cc_mixed,
                k_nue_cc_out_fv,
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

}; // End Class Utility


#endif
