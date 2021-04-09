#ifndef UTILITYPLOTTER_H
#define UTILITYPLOTTER_H

#include "Utility.h"
#include "WienerSVD.h"
#include "SliceContainer.h"
#include "SelectionCuts.h"

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

    std::vector<std::string> var_labels_xsec = {};

    std::vector<std::string> var_labels_events = {};

    std::vector<std::string> var_labels_eff = {};

    std::string smear_hist_name;
    
    std::vector<std::string> vars = {};

    double xsec_scale;

    // enum for histogram vars
    enum TH1D_xsec_var_vars {
        k_var_integrated,     // Total X-Section
        k_var_recoX,          // X-Sec as a function of a Reconstructed variable
        k_var_trueX,          // X-Sec as a function of a True variable
        k_TH1D_xsec_var_MAX
    };

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
    void GetFitResult(double &mean, double &sigma, float bin_lower_edge, float bin_upper_edge, TTree* tree, bool save_hist, bool &converged, bool draw_fit_results, std::string var);
    // -------------------------------------------------------------------------
    // Caller function to optimise the bins
    void OptimiseBins();
    // -------------------------------------------------------------------------
    // Caller function to plot variables by reconstructed shower energy bin
    void PlotVarbyRecoBin();
    // -------------------------------------------------------------------------
    // Function that plots a variable in different bin ranges
    void PlotQuery(float bin_lower_edge, float bin_upper_edge, TTree* tree, std::string xvar, std::string reco_var, std::string true_var);
    // -------------------------------------------------------------------------
    // Get the integrated flux, draw threshold line for technote
    void PlotIntegratedFluxwithThrehold();
    // -------------------------------------------------------------------------
    // Function to plot a number of true variables at the start of the selection
    // e.g. the hit purity and pion mommentum
    void PlotTrueVar();
    // -------------------------------------------------------------------------
    // Function to save a few 2D histograms
    void Save2DHists(const char* printname, TH2D* hist);
    // -------------------------------------------------------------------------
    // Row Normalise 2D histogram
    void RowNorm(TH2D* hist);
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
    // See how a response matrix made with a different model impacts the 
    // smeared MC truth prediction
    void CompareSmearing();
    // -------------------------------------------------------------------------
    // Compare the cross sections for different models in true space
    void CompareUnfoldedModels();
    // -------------------------------------------------------------------------
    // Compare the true cross section of fake data to the extracted cross section
    // in reco space
    void CompareFakeDataReco();
    // -------------------------------------------------------------------------
     // Compare the true cross section of fake data to the extracted cross section
    // in true space
    void CompareFakeDataTrue();
    // -------------------------------------------------------------------------
    // Compare the data to a number of alternative models
    void CompareTotalCrossSec();
    // -------------------------------------------------------------------------
    // Compare the total cross section for each fake data model
    void CompareFakeTotalCrossSec();
    // -------------------------------------------------------------------------
    // Compare the total data cross sections extracted for each model
    void CompareTotalDataCrossSections();
    // -------------------------------------------------------------------------
    // Compare the unfolded data cross section extracted with different models
    void CompareUnfoldedDataCrossSections();
    // -------------------------------------------------------------------------
    // Get the response matrix from file and save it for the technote
    void SaveResponseMatrix();
    // -------------------------------------------------------------------------
    // Check if the pi0 tune is covered by the genie systeamtics
    void CheckPi0Coverage();
    // -------------------------------------------------------------------------
    // Compare the MCC9 result to MCC8
    void CompareMCC8Result();
    // -------------------------------------------------------------------------
    // Compare generator predictions to the data
    void ForwardFoldedGeneratorComparison();
    // -------------------------------------------------------------------------
    // Compare the total cross section generator comparisons
    void CompareGeneratorTotalCrossSec();
    // -------------------------------------------------------------------------
    // Compare unfoldef generators
    void CompareGeneratorUnfoldedModels();
    // -------------------------------------------------------------------------
    // Compare the True generator pi0 distributions
    void CompareGeneratorPi0();
    // -------------------------------------------------------------------------


}; // End Class UtilityPlotter

#endif