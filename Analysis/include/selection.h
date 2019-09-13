#ifndef SELECTION_h
#define SELECTION_h

#include "utility.h"

/* 
Class is the main controller of the selection. Main.cxx will invoke the intilalisation
of this class. After the initialisation, the make_selection() function will be
run.

*/

namespace xsecSelection {

    class selection {

        private:

            // Bools
            bool bool_use_mc{false};
            bool bool_use_ext{false};
            bool bool_use_data{false};
            bool bool_use_dirt{false};
            bool bool_use_variation{false};
            bool slim{false};               // Flag to decide whether to make, fill and plot histograms

            // TFiles
            TFile * f_mc;        // The MC file
            TFile * f_data;      // The data file
            TFile * f_ext;       // The ext file
            TFile * f_dirt;      // The dirt file
            TFile * f_variation; // The variation file

            // TTrees
            TTree * mytree; // MC TPC Object Tree
            TTree * optree; // MC Optical Tree
            TTree * mctree; // MC Tree
            TTree * mctruth_counter_tree; // Tree with counts for MC variables

            TTree * data_tree;   // Data TPC Object Tree
            TTree * data_optree; // Data Optical Tree

            TTree * ext_tree;   // EXT TPC Object Tree
            TTree * ext_optree; // EXT Optical Tree

            TTree * dirt_tree;   // Dirt TPC Object Tree
            TTree * dirt_optree; // Dirt Optical Tree

            // Class Instances
            utilityNS::utility _utility_instance;

            // Variables -------------------------------------------------------

            // Selection Cut Values
            double _x1;
            double _x2;
            double _y1;
            double _y2;
            double _z1;
            double _z2;
            double flash_pe_threshold;
            double flash_time_start;
            double flash_time_end;
            double tolerance;
            double shwr_nue_tolerance;
            double trk_nue_tolerance;
            double shwr_hit_threshold;
            double shwr_hit_threshold_collection;
            double tolerance_open_angle_min;
            double tolerance_open_angle_max;
            double tolerance_dedx_min;
            double tolerance_dedx_max;
            double dist_tolerance;
            double pfp_hits_length_tolerance;
            double ratio_tolerance;
            bool   detector_variations;
            std::vector<double> fv_boundary_v;

            // TPC Object Containers
            std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v       = nullptr;  // MC
            std::vector<xsecAna::TPCObjectContainer> * data_tpc_object_container_v  = nullptr;  // Data
            std::vector<xsecAna::TPCObjectContainer> * ext_tpc_object_container_v   = nullptr;  // EXT
            std::vector<xsecAna::TPCObjectContainer> * dirt_tpc_object_container_v  = nullptr;  // Dirt

            // Optical Containers
            std::vector<std::vector<double>> mc_largest_flash_v_v;    // MC container for largest flashes
            std::vector<std::vector<double>> data_largest_flash_v_v;  // Data container for largest flashes
            std::vector<std::vector<double>> ext_largest_flash_v_v;   // EXT container for largest flashes
            std::vector<std::vector<double>> dirt_largest_flash_v_v;  // Dirt container for largest flashes

            // Counter variables
           int tree_total_entries;      // MC
           int data_tree_total_entries; // Data
           int ext_tree_total_entries;  // EXT
           int dirt_tree_total_entries;  // Dirt

            // MC Truth Counter Tree Variables
            int mc_nue_cc_counter      = 0;
            int mc_nue_nc_counter      = 0;
            int mc_numu_cc_counter     = 0;
            int mc_numu_nc_counter     = 0;
            int mc_nue_cc_counter_bar  = 0;
            int mc_numu_cc_counter_bar = 0;
            int mc_nue_nc_counter_bar  = 0;
            int mc_numu_nc_counter_bar = 0;
            double mc_nu_energy        = 0;
            double mc_nu_momentum      = 0;
            int mc_nu_id               = -1;
            int mc_nu_id_var           = -1;
            double mc_nu_vtx_x         = -999;
            double mc_nu_vtx_y         = -999;
            double mc_nu_vtx_z         = -999;
            double mc_nu_vtx_x_var     = -999;
            double mc_nu_vtx_y_var     = -999;
            double mc_nu_vtx_z_var     = -999;

            double mc_nu_dir_x         = -999;
            double mc_nu_dir_y         = -999;
            double mc_nu_dir_z         = -999;
            double mc_ele_dir_x        = -999;
            double mc_ele_dir_y        = -999;
            double mc_ele_dir_z        = -999;
            double mc_ele_energy       = 0;
            double mc_ele_momentum     = 0;
            bool has_pi0               = false;
            double mc_nu_time          = -1;

            int mc_nu_num_particles         = 0;
            int mc_nu_num_charged_particles = 0;

            // Nue Counter Variables
            int _mc_nue_cc_counter      = 0;
            int _mc_nue_cc_counter_bar  = 0;
            int _mc_numu_cc_counter     = 0;
            int _mc_numu_cc_counter_bar = 0;
            int _mc_nue_nc_counter      = 0;
            int _mc_nue_nc_counter_bar  = 0;
            int _mc_numu_nc_counter     = 0;
            int _mc_numu_nc_counter_bar = 0;
            
            // Vector to check if in TPC
            std::vector<bool> true_in_tpc_v;
            
            int total_mc_entries_inFV_nue         = 0;
            int total_mc_entries_inFV_nue_bar     = 0;
            int total_mc_entries_inFV_numu_cc     = 0;
            int total_mc_entries_inFV_nue_nc      = 0;
            int total_mc_entries_inFV_numu_nc     = 0;
            int total_mc_entries_inFV_numu_cc_bar = 0;
            int total_mc_entries_inFV_nue_nc_bar  = 0;
            int total_mc_entries_inFV_numu_nc_bar = 0;
            int total_mc_entries_inFV             = 0;



        public:
            selection() = default;
            // -----------------------------------------------------------------
            // Allow the utility class to modify the private member data of this class
            friend class utilityNS::utility;
            // -----------------------------------------------------------------
            // Function to Initialise the selection class ----------------------
            void Initialise( const char * mc_file,
                             const char * ext_file,
                             const char * data_file,
                             const char * dirt_file,
                             const char * variation_file,
                             const std::vector<double> _config,
                             bool _slim);
            // -----------------------------------------------------------------
           
            // Main function for selection -------------------------------------
            void make_selection();
            // -----------------------------------------------------------------

    }; // END Class

} // END NAMESPACE xsecSelection
#endif