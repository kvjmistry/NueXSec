#ifndef SYSTEMATICSHELPER_h
#define SYSTEMATICSHELPER_h

#include "utility.h"

// Class for making plots for systematic studies
class SystematicsHelper{

    public:
    // Default constructor
    SystematicsHelper(){};
    
    // The output file
    TFile* f_nuexsec;

    // Class instances
    utility _util;

    // Variables
    int run{0}, subrun{0}, event{0};
    
    TTree * tree;

    std::string run_period;


    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char *run_period, const char * xsec_file_in, utility _utility);
    // -------------------------------------------------------------------------
    // Function to loop over events and calculate the cross section
    void LoopEvents(); 
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    private:

    // Here we create the trees 



}; // End Class SystematicsHelper

#endif