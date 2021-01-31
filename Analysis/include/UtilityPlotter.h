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

    std::string *classification = NULL; // The classification of the event
    
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
    std::vector<unsigned short> *weightsPPFX = NULL ;




    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(Utility _utility);
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
    // Get the integrated flux, draw threshold line for technote
    void PlotIntegratedFluxwithThrehold();
    // -------------------------------------------------------------------------
    // Function to plot a number of true variables at the start of the selection
    // e.g. the hit purity and pion mommentum
    void PlotTrueVar();
    // -------------------------------------------------------------------------
    // Similar function to the slice container classifier, just re-implement it here for easier use
    std::pair<std::string, int> Classify(float true_nu_vtx_sce_x, float true_nu_vtx_sce_y, float true_nu_vtx_sce_z, int nu_pdg, int ccnc, float nu_purity_from_pfp, int npi0);
    // -------------------------------------------------------------------------
    // Function to save a few 2D histograms
    void Save2DHists(const char* printname, TH2D* hist);
    // -------------------------------------------------------------------------
    // Column Normalise 2D histogram
    void ColumnNorm(TH2D* hist);
    // -------------------------------------------------------------------------
    // Study the ppfx weights for each event classification
    void StudyPPFXWeights();
    // -------------------------------------------------------------------------
    // Study the efficeincy for run 1 and run 3
    void CompareEfficiency();
    // -------------------------------------------------------------------------
    // Get the efficiency from the ttree and set the errors
    void PopulateEff(TH1D* h_eff, TH1D *h_pur, TH1D* h_eff_clone, const char* input_file);
    // -------------------------------------------------------------------------
    // Compare the efficiency in the standard det var CV and intrinsic nue det var sample
    void CompareDetVarEfficiency();
    // -------------------------------------------------------------------------
    // Function to compare how making a smearing matrix with a different model impacts
    // the measured cross section.
    void TestModelDependence();
    // -------------------------------------------------------------------------
    // Compare the extracted data cross section for different models
    void CompareDataCrossSections();
    // -------------------------------------------------------------------------



}; // End Class UtilityPlotter

#endif