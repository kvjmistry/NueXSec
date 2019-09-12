#ifndef SELECTION_h
#define SELECTION_h

#include "main.h"
#include "utility.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

namespace xsecSelection {

    class selection {

        private:

            // Bools
            bool bool_use_mc{false};
            bool bool_use_ext{false};
            bool bool_use_data{false};
            bool bool_use_dirt{false};
            bool bool_use_variation{false};

            // TFiles
            TFile * f_mc;        // The MC file
            TFile * f_data;      // The data file
            TFile * f_ext;       // The ext file
            TFile * f_dirt;      // The dirt file
            TFile * f_variation; // The variation file

            // TTrees
            TTree * mytree; // TPC Object Tree
            TTree * optree; // Optical Tree
            TTree * mctree; // MC Tree
            TTree * mctruth_counter_tree; // Tree with counts for MC variables

            // Class Instances
            utilityNS::utility _utility_instance;

        public:
            selection() = default;

            friend class utilityNS::utility;
           
            // Function to Initialise the selection class
            void Initialise( const char * mc_file,
                             const char * ext_file,
                             const char * data_file,
                             const char * dirt_file,
                             const char * variation_file,
                             const std::vector<double> config
                             );
            
           
            // Main function for selection
            void make_selection();

            



    }; // END Class

} // END NAMESPACE xsecSelection
#endif