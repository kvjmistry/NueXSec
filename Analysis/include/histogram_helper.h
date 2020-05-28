#ifndef HISTOGRAM_HELPER_h
#define HISTOGRAM_HELPER_h

#include "selection_cuts.h"
#include "SliceContainer.h"

// Class for filling and saving histograms. 
class histogram_helper{

    public:
    // Default constructor
    histogram_helper(){};
    
    // Destructor 
    ~histogram_helper(); 

    // The output file
    TFile* f_nuexsec;

    // Class instances
    utility _util;
    selection_cuts _scuts;
    int _type{1};

    // weight variable (will equal multiple of all weights)
    double weight{1.0};

    bool weight_tune = true; // Apply genie tune weight
    bool weight_ppfx = true; // Apply ppfx cv weight

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type, const char *run_period, const char * file_out, int weight_cfg );
    // -------------------------------------------------------------------------
    // Initialise histograms
    void InitHistograms();
    // -------------------------------------------------------------------------
    // Function to make the directory hirarchy
    void MakeDirectory();
    // -------------------------------------------------------------------------
    // Function to fill the reco variables
    void FillHists(int type, int classification_index, std::string interaction, int _par_type, int cut_index, SliceContainer SC, double weight);
    // -------------------------------------------------------------------------
    // Function to write the histograms to a file
    void WriteReco(int type);
    // -------------------------------------------------------------------------
    // As above, but for histograms by particle type
    void WriteRecoPar(int type);
    // -------------------------------------------------------------------------
    // Function to fill the true neutrino in FV graphs
    void FillTEfficiency(int cut_index, std::string classification, SliceContainer SC, double weight);
    // -------------------------------------------------------------------------
    // Function to write the TEfficiency Graphs to file
    void WriteTEfficiency();
    // -------------------------------------------------------------------------
    // Function to write the truth histograms to file
    void WriteTrue();
    // -------------------------------------------------------------------------
    // Function to write the flash histograms
    void WriteFlash();
    // -------------------------------------------------------------------------
    // Function to write the interaction histograms
    void WriteInteractions();
    // -------------------------------------------------------------------------
    // Function to write the 2D signal vs Background Histograms
    void Write_2DSigBkgHists();
    // -------------------------------------------------------------------------


    private:

    // Here we create the histograms
   
    // vector of histograms to make, indexed by enums
    std::vector<std::vector<std::vector<TH1D*>>> TH1D_hists; 

    // vector of histograms to make, indexed by enums -- for particle type
    std::vector<std::vector<std::vector<TH1D*>>> TH1D_hists_particle; 

    // Histograms for the efficiency plot
    std::vector<TH1D*> TEfficiency_hists;

    // True histograms
    std::vector<TH1D*> TH1D_true_hists;
    std::vector<TH2D*> TH2D_true_hists;
    
    // Flash Histograms
    std::vector<TH1D*> TH1D_flash_hists;

    // Interaction Histograms
    std::vector<TH1D*> TH1D_interaction_hists;

    // 2D histograms for Signal and Background Rejection
    std::vector<std::vector<TH2D*>> TH2D_hists;
    
    // enum for histogram vars
    enum TH1D_hist_vars {
        k_TH1D_MAX
    };

    // enum for histogram vars
    enum TH1D_par_hist_vars {
        k_reco_dEdx_cali_y_plane_par,                                           // cali dEdx in the collection plane
        k_reco_shr_tkfit_dedx_y_par,                                                // dEdx in the collection plane with trackfit 1x4 cm box
        k_TH1D_par_MAX
    };

    enum TH1D_true_hist_vars {
        k_true_nue_theta,     // True nue in BNB theta coordinates (up from beam dir)
        k_true_nue_phi,       // True nue in BNB phi coordinates (around beam dir)
        k_true_nue_angle,     // True nue angle from numi beamline 
        k_true_nue_px,        // True nue px
        k_true_nue_py,        // True nue py
        k_true_nue_pz,        // True nue pz
        k_true_nue_e,         // True nue energy
        k_true_nue_p,         // True nue momentum
        k_true_vtx_x,         // True Vertex X
        k_true_vtx_y,         // True Vertex Y
        k_true_vtx_z,         // True Vertex Z
        k_true_vtx_x_sce,     // True Vertex X Space Charge Corrected
        k_true_vtx_y_sce,     // True Vertex Y Space Charge Corrected
        k_true_vtx_z_sce,     // True Vertex Z Space Charge Corrected
        k_TH1D_true_MAX
    };

    enum TH1D_flash_hist_vars {
        k_flash_time,
        k_flash_pe,
        k_TH1D_flash_MAX
    };

    enum TH2D_true_hist_vars {
        k_true_nue_phi_theta, 
        k_true_nue_energy_theta,
        k_true_nue_energy_phi,
        k_true_nue_energy_angle,
        k_true_nue_vtx_z_y,
        k_true_nue_vtx_z_y_sce,
        k_true_elec_E_reco_elec_E,
        k_true_nu_vtx_x_reco_nu_vtx_x,
        k_true_nu_vtx_y_reco_nu_vtx_y,
        k_true_nu_vtx_z_reco_nu_vtx_z,
        k_TH2D_true_MAX
    };

    // 2D Histograms for separating signal and background
    enum TH2D_reco_hist_vars {
        k_reco_shr_dEdx_shr_dist,
        k_TH2D_reco_MAX
    };

}; // End Class Histogram Helper 

#endif