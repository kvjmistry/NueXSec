#include "../include/UtilityPlotter.h"

// -----------------------------------------------------------------------------
void UtilityPlotter::Initialise(const char *_run_period, Utility _utility, const char* mode){

    std::cout << "Initalising Utility Plotter ..." << std::endl;
    _util = _utility;

    // Set the run period
    run_period = std::string(_run_period);

    // Standard variation mode
    if (std::string(mode) == "default")  {
        f_nuexsec_mc = TFile::Open( Form("files/tree/nuexsec_tree_mc_run%s.root", _run_period ));
        f_nuexsec    = TFile::Open( Form("files/tree/nuexsec_tree_merged_run%s.root", _run_period ));
    }
    else {
        std::cout << "Error I dont know what mode you have configured..." << mode << std::endl;
        exit(1);
    }
        
}
// -----------------------------------------------------------------------------