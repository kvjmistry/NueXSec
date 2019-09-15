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

// Analysis Includes
#include "Passed_Container.h"
#include "../Modules/LinkDef.h"
#include "../../Modules/TpcObjectContainer.h"
#include "../../Modules/ParticleContainer.h"


namespace utilityNS {

    class utility{

        public:
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
        // Checks in FV, makes no cut, but returns a bool
        bool in_fv(double x, double y, double z, std::vector<double> fv_boundary_v);
        // ---------------------------------------------------------------------
        // Get the largest flash vector from the optical tree
        void GetLargestFlashVector(TTree* optical_tree, double flash_time_start, double flash_time_end, int flash_pe_threshold, 
                                   std::vector<std::vector<double>> &largest_flash_v_v, std::vector<std::vector<int>> &optical_list_pe_v, std::vector<std::vector<double>> &optical_list_flash_time_v);
        // ---------------------------------------------------------------------


    }; // End Class Utility
} // End namespace

#endif
