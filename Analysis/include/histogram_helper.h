#ifndef HISTOGRAM_HELPER_h
#define HISTOGRAM_HELPER_h

#include "selection_cuts.h"

// Class for filling and saving histograms. 
class histogram_helper{

    public:
    // Default constructor
    histogram_helper(){};
    
    // Destructor 
    ~histogram_helper(); 

    // The output file
    TFile* f_nuexsec = new TFile("nuexsec.root", "UPDATE");

    // Class instances
    utilityNS::utility _utility_instance;
    selection_cuts _selection_cuts_instance;

    // enum to switch file type 
    enum type {k_mc, k_data, k_ext, k_dirt, k_variation, k_type_MAX}; 

    // For creating histogram names
    std::vector<std::string> type_prefix = {"MC", "Data", "EXT", "Dirt"};

    // enum for vertex zy histograms
    enum vertex {k_vtx_signal, k_vtx_data, k_vtx_ext, k_vtx_mc_ext, k_vertex_MAX};
    std::vector<std::string> vertex_strings = {"Signal", "Data", "EXT", "MC_plus_EXT"};
    
    // enums for cut dirs
    enum enum_cut_dirs {
                k_in_fv,                           // Fiducial volume
                k_cuts_MAX
                };
    
    
    // enums for flash times
    enum flash_time{ k_flash_mc, k_flash_data, k_flash_ext, k_flash_dirt, k_flash_MAX };

    // Cut directory names
    std::vector<std::string> cut_dirs = {
                "In_FV"                                    // Fiducial volume
                };

    
    // Names of the plot types
    std::vector<std::string> plot_types = {
                "Truth",
                "Reco",
                "Optical",
                "Stack"
                };
    
    // Names of the classifications
    std::vector<std::string> classification_dirs = {
                "nue_cc",
                "nue_cc_mixed",
                "nue_cc_out_fv",
                "cosmic",
                "numu_cc",
                "nc",
                "nc_pi0",
                "nc_mixed",
                "unmatched",
                "ext",
                "data",
                "dirt"
                };

    
    enum legend {
                k_nue_cc,
                k_nue_cc_mixed,
                k_nue_cc_out_fv,
                k_cosmic,
                k_numu_cc,
                k_nc,
                k_nc_pi0,
                k_nc_mixed,
                k_unmatched,
                k_leg_ext,
                k_leg_data,
                k_leg_dirt,
                k_classifications_MAX
                };

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise();
    // -------------------------------------------------------------------------
    // Initialise histograms
    void InitHistograms();
    // -------------------------------------------------------------------------
    // Function to make the directory hirarchy
    void MakeDirectory(std::string type);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    private:

    // Here we create the histograms
    // --------------------------- True plots ----------------------------------
    TH1D * h_nue_true_theta  = new TH1D ("h_nue_true_theta","True Nue; Theta [degrees]", 14,    0, 180);
    TH1D * h_nue_true_phi    = new TH1D ("h_nue_true_phi"  ,"True Nue; Phi [degrees]",   14, -180, 180);
    
    TH2D * h_nue_true_theta_phi    = new TH2D("h_nue_true_theta_phi",    "True Nue; Phi [degrees]; Theta [degrees]", 14, -20,100,  14,   20, 140);
    TH2D * h_nue_true_energy_theta = new TH2D("h_nue_true_energy_theta", "True Nue; Energy [GeV]; Theta[degrees]",   20,   0, 10,  14,    0, 180);
    TH2D * h_nue_true_energy_phi   = new TH2D("h_nue_true_energy_phi",   "True Nue; Theta [degrees]; Phi [degrees]", 20,   0, 10,  14, -180, 180);
    TH2D * h_ele_true_energy_theta = new TH2D("h_ele_true_energy_theta", "True e; Theta [degrees]; Theta [degrees]", 20,   0, 10,  14,    0, 180);
    TH2D * h_ele_true_energy_phi   = new TH2D("h_ele_true_energy_phi",   "True e; Theta [degrees]; Phi [degrees]",   20,   0, 10,  14, -180, 180);

    // Optical Plots
    std::vector<TH1D*> h_flash_time_v;
    std::vector<std::vector<TH1D*>> h_largest_flash_z; // Largest flash Z --
   
    // RecoPlots
    std::vector<std::vector<TH1D*>> h_reco_vtx_x; // Reco Vertex X
    std::vector<std::vector<TH1D*>> h_reco_vtx_y; // Reco Vertex Y
    std::vector<std::vector<TH1D*>> h_reco_vtx_z; // Reco Vertex Z
    std::vector<TH2D*> h_reco_vtx_zy;             // Reco Vertex ZY Plane

    // Leading Shower Momentum
    std::vector<std::vector<TH1D*>> h_reco_leading_mom;

    // 2D distance largest flash to reco nu vertex
    std::vector<std::vector<TH1D*>> h_reco_flash_to_vtx_dist;

    // 2D distance shower vertex to reco nu vertex
    std::vector<std::vector<TH1D*>> h_reco_shower_to_vtx_dist;

    // 2D distance track vertex to reco nu vertex
    std::vector<std::vector<TH1D*>> h_reco_track_to_vtx_dist;

    // Leading Shower hits in all planes
    std::vector<std::vector<TH1D*>> h_reco_leading_shower_hits_all_planes;

    // Leading Shower hits in collection
    std::vector<std::vector<TH1D*>> h_reco_leading_shower_hits_collection_plane;

    // Leading Shower opening angle
    std::vector<std::vector<TH1D*>> h_reco_leading_shower_open_angle;

    // dEdx
    std::vector<std::vector<TH1D*>> h_reco_dEdx;

    // Secondary shower to vertex distance (for events with more than 1 shower)
    std::vector<std::vector<TH1D*>> h_reco_secondary_shower_to_vtx_dist;

    // Leading Shower hits per length
    std::vector<std::vector<TH1D*>> h_reco_leading_shower_hits_per_length;

    // Longest track to leading shower length
    std::vector<std::vector<TH1D*>> h_reco_longest_track_leading_shower_length;

    // Track Containment
    std::vector<std::vector<TH1D*>> h_reco_track_contained;

    // Leading shower phi
    std::vector<std::vector<TH1D*>> h_reco_leading_shower_phi;

    // Leading shower theta
    std::vector<std::vector<TH1D*>> h_reco_leading_shower_theta;

    // Leading shower cos theta
    std::vector<std::vector<TH1D*>> h_reco_leading_shower_cos_theta;

    // Leading shower multiplicity
    std::vector<std::vector<TH1D*>> h_reco_shower_multiplicity;

    // Leading track multiplicity
    std::vector<std::vector<TH1D*>> h_reco_track_multiplicity;

}; // End Class Histogram Helper 

#endif