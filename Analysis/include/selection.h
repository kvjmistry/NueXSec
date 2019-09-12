#ifndef SELECTION_h
#define SELECTION_h

#include "main.h"

namespace xsecSelection {

    class selection {

        private:

        // TFiles
        // TFile * f_mc;        // The MC file
        // TFile * f_data;      // The data file
        // TFile * f_ext;       // The ext file
        // TFile * f_dirt;      // The dirt file
        // TFile * f_variation; // The variation file


        public:
            selection() = default;
           
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