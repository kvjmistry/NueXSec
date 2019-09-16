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

    // enum to switch file type 
    enum type {k_mc, k_data, k_ext, k_dirt, k_variation, k_type_MAX}; 

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
    void FillOptical(std::vector<std::vector<double>> optical_list_flash_time_v, int type);
    // -------------------------------------------------------------------------
    // Write Optical information to file
    void WriteOptical(int type);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    private:

    // The output file
    TFile* f_nuexsec = new TFile("nuexsec.root", "UPDATE");

    // Cut directory names
    std::vector<std::string> cut_dirs = {
                "Flash_PE_and_In_Time",                     // Flash PE and In time flash
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
                "Optical",
                "Stack"
                };
    

    // Names of the classifications
    std::vector<std::string> classification_dirs = {
                "nue_cc_mixed",
                "numu_cc_mixed",
                "other_mixed",
                "cosmic",
                "nue_cc_out_fv",
                "nue_cc_out_fv",
                "nue_cc_qe",
                "nue_cc_res",
                "nue_cc_dis",
                "nue_cc_coh",
                "nue_cc_mec",
                "nue_bar_cc_qe",
                "nue_bar_cc_res",
                "nue_bar_cc_dis",
                "nue_bar_cc_coh",
                "nue_bar_cc_mec",
                "numu_cc_qe",
                "numu_cc_res",
                "numu_cc_dis",
                "numu_cc_coh",
                "numu_cc_mec",
                "nc",
                "nc_pi0",
                "unmatched"
                };

    // Class instances
    utilityNS::utility _utility_instance;


    // Here we create the histograms

    // For creating histogram names
    std::vector<std::string> type_prefix = {"MC", "Data", "EXT", "Dirt"};
    
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
    enum flash_time{ k_flash_mc, k_flash_data, k_flash_ext, k_flash_dirt, k_flash_MAX };

}; // End Class Histogram Helper

#endif
