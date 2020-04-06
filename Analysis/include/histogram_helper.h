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
    void FillHists(int type, int classification_index, std::string interaction, int cut_index, SliceContainer SC, double weight);
    // -------------------------------------------------------------------------
    // Function to write the histograms to a file
    void WriteReco(int type);
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

    private:

    // Here we create the histograms
   
    // vector of histograms to make, indexed by enums
    std::vector<std::vector<std::vector<TH1D*>>> TH1D_hists; 

    // Histograms for the efficiency plot
    std::vector<TH1D*> TEfficiency_hists;

    // True histograms
    std::vector<TH1D*> TH1D_true_hists;
    std::vector<TH2D*> TH2D_true_hists;
    
    // Flash Histograms
    std::vector<TH1D*> TH1D_flash_hists;

    // Interaction Histograms
    std::vector<TH1D*> TH1D_interaction_hists;

    // enum for histogram vars
    enum TH1D_hist_vars {
        k_reco_vtx_x,                                                           // Reco Vertex X
        k_reco_vtx_y,                                                           // Reco Vertex Y
        k_reco_vtx_z,                                                           // Reco Vertex Z
        k_reco_vtx_x_sce,                                                       // Reco Vertex X Space Charge Corrected
        k_reco_vtx_y_sce,                                                       // Reco Vertex Y Space Charge Corrected
        k_reco_vtx_z_sce,                                                       // Reco Vertex Z Space Charge Corrected
        k_reco_dEdx_y_plane,                                                    // dEdx on Collection Plane (uncalibrated)
        k_reco_dEdx_cali_y_plane,                                               // dEdx Cali on Collection Plane
        k_reco_leading_mom,                                                     // Leading Shower Momentum
        k_reco_shower_to_vtx_dist,                                              // 2D distance shower vertex to reco nu vertex
        k_reco_track_to_vtx_dist,                                               // 2D distance track vertex to reco nu vertex
        k_reco_leading_shower_hits_all_planes,                                  // Leading Shower hits in all planes
        k_reco_leading_shower_hits_collection_plane,                            // Leading Shower hits in collection
        k_reco_leading_shower_open_angle,                                       // Leading Shower opening angle
        k_reco_secondary_shower_to_vtx_dist,                                    // Secondary shower to vertex distance (for events with more than 1 shower)
        k_reco_leading_shower_hits_per_length,                                  // Leading Shower hits per length
        k_reco_longest_track_leading_shower_length,                             // Longest track to leading shower length
        k_reco_n_track_contained,                                               // Number of Tracks Contained
        k_reco_n_shower_contained,                                              // Number of Showers Contained
        k_reco_leading_shower_phi,                                              // Leading shower phi
        k_reco_leading_shower_theta,                                            // Leading shower theta
        k_reco_leading_shower_cos_theta,                                        // Leading shower cos theta
        k_reco_shower_multiplicity,                                             // Leading shower multiplicity
        k_reco_track_multiplicity,                                              // Leading track multiplicity
        k_reco_topological_score,                                               // Pandora Topological Score
        k_reco_track_shower_dist,                                               // Track shower dist
        k_reco_track_shower_angle,                                              // Track shower angle
        k_reco_hits_ratio,                                                      // Ratio hits from showers to slice
        k_reco_shower_score,                                                    // Shower score
        k_reco_track_score,                                                     // Track score
        k_reco_shower_energy_tot_cali,                                          // Calibrated energy of all the showers
        k_reco_shower_hits,                                                     // Total number of hits for the leading shower
        k_reco_shower_hits_y_plane,                                             // Total number of hits for the leading shower in the collection plane
        k_reco_shr_trkfit_2cm_dEdx,                                             // dE/dx of the leading shower on the Y plane with the track fitting, use first 2 cm
        k_reco_shr_trkfit_2cm_dEdx_y,                                           // dE/dx of the leading shower on the Y plane with the track fitting, use first 2 cm y plane
        k_reco_shr_trkfit_gap05_dEdx,                                           // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 10mm
        k_reco_shr_trkfit_gap05_dEdx_y,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 5mm y plane
        k_reco_shr_trkfit_gap10_dEdx,                                           // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 5mm
        k_reco_shr_trkfit_gap10_dEdx_y,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 10mm y plane
        k_reco_opfilter_beam,                                                   // Common optical filter beam
        k_reco_opfilter_veto,                                                   // Common Optical filter michel veto
        k_reco_softwaretrig,                                                    // Software Trigger
        k_reco_nslice,                                                          // Pandora Slice ID
        k_reco_slclustfrac,                                                     // Reco Fraction of hits in the slice that are fully reconstructed to 3D particles.
        k_reco_cosmicIP,                                                        // Reco Closest distance between shower start and space points associated to tracks flagged as cosmics.
        k_reco_shr_tkfit_dedx_Y,                                                // The dEdx using the trackfit variable
        k_TH1D_MAX
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
        k_flash_time_sid1,   // Slice Id 1 (neutrino candidates)
        k_flash_pe_sid1,
        k_flash_time_sid0,  // Slice Id 0 (non neutrino candidates)
        k_flash_pe_sid0,
        k_TH1D_flash_MAX
    };

    enum TH2D_true_hist_vars {
        k_true_nue_theta_phi, 
        k_true_nue_energy_theta,
        k_true_nue_energy_phi,
        k_true_nue_energy_angle,
        k_true_nue_vtx_z_y,
        k_true_nue_vtx_z_y_sce,
        k_TH2D_true_MAX
    };

    


}; // End Class Histogram Helper 

#endif