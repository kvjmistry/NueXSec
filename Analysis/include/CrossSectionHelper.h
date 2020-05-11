#ifndef CROSSSECTIONHELPER_h
#define CROSSSECTIONHELPER_h

#include "utility.h"

// Class for calculating the cross section and doing systematics
class CrossSectionHelper{

    public:
    // Default constructor
    CrossSectionHelper(){};
    
    // The output file
    TFile* f_nuexsec;

    // Class instances
    utility _util;

     // Scale factors (everything is scaled to data)
    double mc_scale_factor     = 1.0;
    double intime_scale_factor = 1.0;
    double dirt_scale_factor   = 1.0;

    // Variables
    int run{0}, subrun{0}, event{0};
    std::string *classifcation = NULL; // The classification of the event
    
    // Is the event a true signal event in the FV that was not selected?
    // We still need these for the efficiency
    bool gen{false};           
    
    double weight{0.0};        // This is not going to be integer if we already weight the CV

    double true_energy{0.0}, reco_energy{0.0};

    TTree * tree;

    std::string run_period;

    // Cross section Variables -- So far copied from coltons analysis, so these numbers need updating
    double lar_density_mc   = 1.3954;   // Density of Argon in the simulation g/cm3
    double lar_density_data = 1.3836;   // Density of Argon in the simulation g/cm3
    double volume           = 4.1622e7; // Fiducial volume cm3
    double NA               = 6.022140857e23; // Advogadro's number molecule/mol
    double N_nuc            = 40.0;     // Number of argon nuclons
    double m_mol            = 39.95;    // Molar Mass of Argon g/mol

    double N_target_MC{0.0};   // Set in code, total number of targets in MC
    double N_target_Data{0.0}; // Set in code, total number of targets in Data

    // Fluxes need to be defined by reading in the flux file and integrating
    double integrated_flux{0.0};
    double flux_scale_factor{1.0e-4}; // m2 to cm2
    double mc_flux_scale_factor{1.0};
    double data_flux_scale_factor{1.0};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char *run_period, const char * xsec_file_in, utility _utility);
    // -------------------------------------------------------------------------
    // Function to loop over events and calculate the cross section
    void LoopEvents(); 
    // -------------------------------------------------------------------------
    // Function to calculate the cross section
    double CalcCrossSec(double sel, double gen, double sig, double bkg, double flux, double ext, double dirt, double targ);
    // -------------------------------------------------------------------------
    // Function to get the integrated flux
    double GetIntegratedFlux();
    // -------------------------------------------------------------------------
    // Function to get the POT from the flux file
    double GetPOT(TFile* f, bool disp);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    private:

    // Here we create the trees 



}; // End Class CrossSectionHelper

#endif