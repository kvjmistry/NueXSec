#ifndef SELECTION_h
#define SELECTION_h

#include "utility.h"
#include "selection_cuts.h"
#include "histogram_helper.h"
#include "histogram_plotter.h"
#include "SliceContainer.h"
#include "Passed_Container.h"

#include <omp.h>
#include "TThread.h"

/* 
Class is the main controller of the selection. Main.cxx will invoke the intilalisation
of this class. After the initialisation, the make_selection() function will be
run.

*/

namespace xsecSelection {

    class selection {

        private:
    
            // Parallel processing variables
            int nthreads = 4;

            int max_events{-1};

            // POT and triggers
            double MC_POT    = 1.0;    // The MC POT
            double Dirt_POT  = 1.0;    // The dirt POT
            double Data_POT  = 1.0;    // The Data POT (tortgt_wcut)
            double Data_trig = 1.0;    // The data HW Triggers (EA9CNT_wcut)
            double EXT_trig  = 1.0;    // The Number of EXT HW Triggers ( EXT_NUMIwin_FEMBeamTriggerAlgo )


            // Scale factors (veverything is scaled to data)
            double mc_scale_factor     = 1.0;
            double intime_scale_factor = 1.0;
            double dirt_scale_factor   = 1.0;

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
            TTree * mc_tree;      // MC   Tree
            TTree * data_tree;    // Data Tree
            TTree * ext_tree;     // EXT  Tree
            TTree * dirt_tree;    // Dirt Tree

            // Class Instances
            utility _util;
            selection_cuts    _scuts;
            std::vector<histogram_helper>  _hhelper; // One for each type e.g MC, Data, EXT..

            // Variables -------------------------------------------------------

            // Slice Containers
            SliceContainer mc_SC;
            SliceContainer data_SC;
            SliceContainer ext_SC;
            SliceContainer dirt_SC;

            // Selection Cut Values
            bool   detector_variations;

            // Vector containing the cut names
            std::vector<std::string> cut_names = {
                "In FV"                                    // Fiducial volume
            };

            // Passed Containers
            std::vector<Passed_Container> mc_passed_v;   // MC Passed Container
            std::vector<Passed_Container> data_passed_v; // Data Passed Container
            std::vector<Passed_Container> ext_passed_v;  // EXT Passed Container
            std::vector<Passed_Container> dirt_passed_v; // Dirt Passed Container
            
            // Counter Containers
            std::vector<std::vector<double>> counter_v;

            // Counter variables
            int mc_tree_total_entries{0};    // MC
            int data_tree_total_entries{0};  // Data
            int ext_tree_total_entries{0};   // EXT
            int dirt_tree_total_entries{0};  // Dirt

            // Total nue's in the cryostat
            int tot_true_cryo_nues;

        public:
            // -----------------------------------------------------------------
            // Allow the utility class to modify the private member data of this class
            friend class utility;
            // -----------------------------------------------------------------
            // Function to Initialise the selection class ----------------------
            void Initialise( const char * mc_file,
                             const char * ext_file,
                             const char * data_file,
                             const char * dirt_file,
                             const char * variation_file,
                             const std::vector<double> _config,
                             bool _slim,
                             int num_events,
                             const char * run_period);
            // -----------------------------------------------------------------
            // Main function for selection
            void make_selection();
            // -----------------------------------------------------------------
            // Template code to apply selection cuts
            bool ApplyCuts(int type, int ievent,std::vector<std::vector<double>> &counter_v,
                           std::vector<Passed_Container> &passed_v, SliceContainer &SC);
            // -----------------------------------------------------------------
            // Function to save all written histograms to file
            void SavetoFile();
            // -----------------------------------------------------------------
            // Make final histogram plots
            void MakeHistograms(const char * hist_file_name, const char *run_period);
            // -----------------------------------------------------------------

    }; // END Class

} // END NAMESPACE xsecSelection
#endif