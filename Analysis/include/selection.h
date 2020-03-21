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

            // Decide whether to make a run subrun event filelist for the selected data events
            bool make_list{true};

            // TFiles
            TFile * f_mc;        // The MC file
            TFile * f_data;      // The data file
            TFile * f_ext;       // The ext file
            TFile * f_dirt;      // The dirt file
            TFile * f_variation; // The variation file
            TFile * f_flux_weights; // The file with the flux weights

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

            int _run_period;

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
            void MakeHistograms(const char * hist_file_name, const char *run_period, const std::vector<double> _config);
            // -----------------------------------------------------------------

    }; // END Class

} // END NAMESPACE xsecSelection
#endif