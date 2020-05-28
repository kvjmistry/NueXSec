#ifndef SELECTION_h
#define SELECTION_h

#include "utility.h"
#include "selection_cuts.h"
#include "histogram_helper.h"
#include "histogram_plotter.h"
#include "SliceContainer.h"
#include "Passed_Container.h"
#include "TreeHelper.h"

#include <omp.h>
#include "TThread.h"

/* 
Class is the main controller of the selection. Main.cxx will invoke the intilalisation
of this class. After the initialisation, the MakeSelection() function will be
run.

*/

namespace xsecSelection {

    class selection {

        private:
    
            // Parallel processing variables
            int nthreads = 4;

            int verbose = 1;

            int max_events{-1};

            int _weight_cfg{1};

            // Scale factors (everything is scaled to data)
            double mc_scale_factor     = 1.0;
            double intime_scale_factor = 1.0;
            double dirt_scale_factor   = 1.0;

            // Bools
            bool bool_use_mc{false};
            bool bool_use_ext{false};
            bool bool_use_data{false};
            bool bool_use_dirt{false};
            bool slim{false};               // Flag to decide whether to make, fill and plot histograms

            // Decide whether to make a run subrun event filelist for the selected data events
            bool make_list{true};

            // TFiles
            TFile * f_mc;        // The MC file
            TFile * f_data;      // The data file
            TFile * f_ext;       // The ext file
            TFile * f_dirt;      // The dirt file
            TFile * f_flux_weights; // The file with the flux weights

            // TTrees
            TTree * mc_tree;      // MC   Tree
            TTree * mc_truth_tree; // MC Truth Tree
            TTree * data_tree;    // Data Tree
            TTree * ext_tree;     // EXT  Tree
            TTree * dirt_tree;    // Dirt Tree

            TTree * eff_tree_out; // TTree to write to the efficiency vec to file for plotting

            // Class Instances
            utility _util;
            selection_cuts    _scuts;
            std::vector<histogram_helper>  _hhelper; // One for each type e.g MC, Data, EXT..
            std::vector<TreeHelper>        _thelper; // One for each type e.g MC, Data, EXT..

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
            
            // Counter Container
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
                             const char * mc_file_out,
                             const char * ext_file_out,
                             const char * data_file_out,
                             const char * dirt_file_out,
                             const char * mc_tree_file_name_out,
                             utility _utility,
                             bool _slim,
                             int num_events,
                             const char * run_period,
                             int _verbose,
                             int weight_cfg);
            // -----------------------------------------------------------------
            // Main function for selection
            void MakeSelection();
            // -----------------------------------------------------------------
            // Template code to apply selection cuts
            bool ApplyCuts(int type, int ievent,std::vector<std::vector<double>> &counter_v,
                           std::vector<Passed_Container> &passed_v, SliceContainer &SC);
            // -----------------------------------------------------------------
            // Function to save all written histograms to file
            void SavetoFile();
            // -----------------------------------------------------------------
            // Function to implement counters and fillings of histograms for each cut
            void SelectionFill(int type, SliceContainer &SC, std::pair<std::string, int> classification, std::string interaction, std::pair<std::string, int> par_type, int cut_index, std::vector<std::vector<double>> &counter_v);
            // -----------------------------------------------------------------

    }; // END Class

} // END NAMESPACE xsecSelection
#endif