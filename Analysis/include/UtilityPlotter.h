#ifndef UTILITYPLOTTER_H
#define UTILITYPLOTTER_H

#include "Utility.h"

// Class for making plots of generic things, such as run vs run comparisons and 
// separate stuies that people want me to do for the analsis


class UtilityPlotter{

    public:
    // Default constructor
    UtilityPlotter(){};
    
    // The output file
    std::vector<TFile*> f_vars;

    // Input file(s)
    TFile *f_nuexsec_mc; // just mc file
    TFile *f_nuexsec;    // merged file

    // Class instances
    Utility _util;

    // Variables
    TTree * tree;
    std::string run_period;

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char *run_period, Utility _utility, const char* mode);
    // -------------------------------------------------------------------------

}; // End Class UtilityPlotter

#endif