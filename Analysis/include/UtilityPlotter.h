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
    TFile *f_nuexsec;    // merged file

    // Class instances
    Utility _util;

    // Variables
    TTree * tree;
    std::string run_period;

    std::string *classifcation = NULL; // The classification of the event
    
    // Is the event a true signal event in the FV that was not selected?
    // We still need these for the efficiency
    bool gen{false};           
    
    double weight{0.0};        // This is not going to be integer if we already weight the CV

    double true_energy{0.0}, reco_energy{0.0};
    float shr_energy_cali{0.0};
    float elec_e{0.0};
    float ppfx_cv{1.0};
    float weightSplineTimesTune{1.0};
    float numi_ang{0.0};
    int nu_pdg{0};
    int shr_bkt_pdg{0};
    float shr_bkt_purity{0.0};
    float shr_bkt_completeness{0.0};
    float shr_bkt_E{0.0}; // energy of truth matched particle to the leading shower
    std::vector<float> *all_shr_hits = NULL;
    std::vector<float> *all_shr_energies = NULL;




    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char *run_period, Utility _utility, const char* mode);
    // -------------------------------------------------------------------------
    // Initialsise the TTree
    void InitTree();
    // -------------------------------------------------------------------------
    // Function to compare the number of hits of the leading shower to see if it
    // really is the shower with the most energy
    void CompareHitstoEnergy();
    // -------------------------------------------------------------------------
    // Function to see how many of the signal events we dont select the electron
    void CompareSignalPurity();
    // -------------------------------------------------------------------------
    // Function to optimise binning for the selection
    void GetFitResult(double &mean, double &sigma, float bin_lower_edge, float bin_upper_edge, TTree* tree, bool save_hist, bool &converged, bool draw_fit_results);
    // -------------------------------------------------------------------------
    // Caller function to optimise the bins
    void OptimiseBins();
    // -------------------------------------------------------------------------
    // Caller function to plot variables by reconstructed shower energy bin
    void PlotVarbyRecoBin();
    // -------------------------------------------------------------------------
    // Function that plots a variable in different bin ranges
    void PlotQuery(float bin_lower_edge, float bin_upper_edge, TTree* tree, std::string variable);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------



}; // End Class UtilityPlotter

#endif