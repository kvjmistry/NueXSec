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
                k_pandora_output,                  // The output TPC Objects from Pandora
                k_in_time,                         // In time flash
                k_flash_pe,                        // Flash PE and In time flash
                k_has_nue,                         // Pandora Nue 
                k_in_fv,                           // Fiducial volume
                k_vtx_to_flash,                    // Vertex to flash
                k_shwr_nue_dist,                   // Distance between pfp shower and nue object
                k_trk_nue_dist,                    // Distance between pfp track and nue object
                k_shwr_hit_threshold,              // Hit threshold for at least one shower
                k_shwr_hit_threshold_collection,   // Hit threshold for at least one shower on collection plane
                k_shwr_open_angle,                 // Tolerance for leading shower open angle
                k_shwr_dedx,                       // Tolerance for dedx of leading shower
                k_dist_nue_vtx,                    // Tolerance for distance from the reco nue vtx for TPCO w/ >3 showers
                k_pfp_hits_length,                 // Tolerance for hits/length
                k_longest_trk_leading_shwr_length, // Tolerance for longest track length / leading shower length
                k_trk_contained,                   // Contained Tracks
                k_cuts_MAX
                };
    
    
    // enums for flash times
    enum flash_time{ k_flash_mc, k_flash_data, k_flash_ext, k_flash_dirt, k_flash_MAX };

    // Cut directory names
    std::vector<std::string> cut_dirs = {
                "Pandora_Output",                           // The output TPC Objects from Pandora
                "In_Time",                                  // In time flash
                "Flash_PE",                                 // Flash PE 
                "Has_Pandora_Nue",                          // Pandora Nue 
                "In_FV",                                    // Fiducial volume
                "Vertex_to_Flash",                          // Vertex to flash
                "Shower_to_Nue_Dist",                       // Distance between pfp shower and nue object
                "Track_to_Nue_Dist",                        // Distance between pfp track and nue object
                "Shower_Hit_Threshold",                     // Hit threshold for at least one shower
                "Shower_Hit_Threshold_Collection_Plane",    // Hit threshold for at least one shower on collection plane
                "Shower_Open_Angle",                        // Tolerance for leading shower open angle
                "Shower_dEdx",                              // Tolerance for dedx of leading shower
                "Dist_Nue_Vtx_Secondary_Showers",           // Tolerance for distance from the reco nue vtx for TPCO w/ >3 showers
                "Hits_per_Length",                          // Tolerance for hits/length
                "Longest_Track_over_Leading_Shower_Length", // Tolerance for longest track length / leading shower length
                "Track_Is_Contained"                        // Contained Tracks
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
    // Fills histograms for truth variables
    void FillMCTruth( double mc_nu_energy,  double mc_nu_momentum,  int mc_nu_id, bool in_tpc,
                      double mc_nu_vtx_x,   double mc_nu_vtx_y,  double mc_nu_vtx_z,
                      double mc_nu_dir_x,   double mc_nu_dir_y,  double mc_nu_dir_z,
                      double mc_ele_dir_x,  double mc_ele_dir_y, double mc_ele_dir_z,
                      double mc_ele_energy, double mc_ele_momentum );
    // -------------------------------------------------------------------------
    // Write MC Truth Variables to file
    void WriteMCTruth(std::string type);
    // -------------------------------------------------------------------------
    // Fill the optical information histograms
    void FillOptical(std::vector<std::vector<double>> &optical_list_flash_time_v, int type);
    // -------------------------------------------------------------------------
    // Write Optical information to file
    void WriteOptical(int type);
    // -------------------------------------------------------------------------
    // Function to return the enum index of the classification type. This will be used
    // to fill the appropriate histogram
    int IndexOfClassification(std::string tpco_id);
    // -------------------------------------------------------------------------
    // Fill reco variables e.g vertex x,y,z, dEdx
    void FillReco(int classification_index, int cut_index, const xsecAna::TPCObjectContainer &tpc_obj, int leading_shower_index, std::vector<double> &largest_flash_v, std::vector<double> fv_boundary_v, std::string type);
    // -------------------------------------------------------------------------
    // Fill reco variables e.g vertex x,y,z, dEdx
    void WriteReco(int type);
    // -------------------------------------------------------------------------
    // Fill reco vertex plot for ZY plane
    void FillRecoVtxZY(int classification_index, const xsecAna::TPCObjectContainer &tpc_obj);
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
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