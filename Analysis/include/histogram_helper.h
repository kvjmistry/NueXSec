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

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type, const char *run_period );
    // -------------------------------------------------------------------------
    // Initialise histograms
    void InitHistograms();
    // -------------------------------------------------------------------------
    // Function to make the directory hirarchy
    void MakeDirectory();
    // -------------------------------------------------------------------------
    // Function to fill the reco variables
    void FillReco(int type, int classification_index, int cut_index, SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Function to write the histograms to a file
    void WriteReco(int type);
    // -------------------------------------------------------------------------
    // Function to fill the true neutrino in FV graphs
    void FillTEfficiency(int cut_index, std::string classification, SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Function to write the TEfficiency Graphs to file
    void WriteTEfficiency();
    // -------------------------------------------------------------------------

    private:

    // Here we create the histograms
    // --------------------------- True plots ----------------------------------
    // TH1D * h_nue_true_theta  = new TH1D ("h_nue_true_theta","True Nue; Theta [degrees]", 14,    0, 180);
    // TH1D * h_nue_true_phi    = new TH1D ("h_nue_true_phi"  ,"True Nue; Phi [degrees]",   14, -180, 180);
    
    // TH2D * h_nue_true_theta_phi    = new TH2D("h_nue_true_theta_phi",    "True Nue; Phi [degrees]; Theta [degrees]", 14, -20,100,  14,   20, 140);
    // TH2D * h_nue_true_energy_theta = new TH2D("h_nue_true_energy_theta", "True Nue; Energy [GeV]; Theta[degrees]",   20,   0, 10,  14,    0, 180);
    // TH2D * h_nue_true_energy_phi   = new TH2D("h_nue_true_energy_phi",   "True Nue; Theta [degrees]; Phi [degrees]", 20,   0, 10,  14, -180, 180);
    // TH2D * h_ele_true_energy_theta = new TH2D("h_ele_true_energy_theta", "True e; Theta [degrees]; Theta [degrees]", 20,   0, 10,  14,    0, 180);
    // TH2D * h_ele_true_energy_phi   = new TH2D("h_ele_true_energy_phi",   "True e; Theta [degrees]; Phi [degrees]",   20,   0, 10,  14, -180, 180);
   
    // vector of histograms to make, indexed by enums
    std::vector<std::vector<std::vector<TH1D*>>> TH1D_hists; 

    // Histograms for the efficiency plot
    std::vector<TH1D*> TEfficiency_hists;
    TH1D* h_true_nu_E; // The true neutrino energy for nu's in the cryostat

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
        k_reco_track_contained,                                                 // Track Containment
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
        k_TH1D_MAX
    };


}; // End Class Histogram Helper 

#endif