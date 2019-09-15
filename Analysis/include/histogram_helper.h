#ifndef HISTOGRAM_HELPER_h
#define HISTOGRAM_HELPER_h

#include "selection_cuts.h"

// Class for filling and saving histograms. 
class histogram_helper{

    private:

    // enum to switch file type 
    enum type {k_mc, k_data, k_ext, k_dirt, k_variation}; 

    // The output file
    TFile* f_nuexsec = new TFile("nuexsec.root", "UPDATE");

    // Cut directory names
    std::vector<std::string> cut_dirs = {
                "Flash_PE_and_In_Time",                     // Flash PE and In time flash
                "Has_Pandora_Nue",                          // Pandora Nue 
                "In_FV",                                    // Fiducial volume
                "Vertex_to_Flash",                          // Vertex to flash
                "Shower_to_Nue_Dist",                       // Distance between pfp shower and nue object
                "Track_to_Nue_Dist",                         // Distance between pfp track and nue object
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

    public:
    // Default constructor
    histogram_helper(){};
    
    // Destructor 
    ~histogram_helper(); 

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise();
    // -------------------------------------------------------------------------
    // Function to make the directory hirarchy
    void MakeDirectory(std::string type);
    // -------------------------------------------------------------------------
    // Fills histograms for truth variables
    void MCTruthInit();
    // -------------------------------------------------------------------------
    

}; // End Class Histogram Helper

#endif
