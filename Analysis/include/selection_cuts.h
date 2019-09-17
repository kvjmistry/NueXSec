#ifndef SELECTION_CUTS_h
#define SELECTION_CUTS_h

#include "utility.h"

// Class for applying the selection cuts. Cut classes will inherit
// from this one and override function
class selection_cuts{

    private:

    public:
    selection_cuts(){};                 // Default constructor

    // Flash variables
    double largest_flash_y{0};
    double largest_flash_z{0};
    double largest_flash_time{0};
    double largest_flash_pe{0};

    // TPC Object Variables
    bool true_in_tpc; // MC Truth -- see if vtx is in tpc
    std::pair<std::string, int> tpc_classification;  // Classification (classification, leading index)

    double tpc_obj_vtx_x = -999.; // TPC Object Vtx X,Y,Z
    double tpc_obj_vtx_y = -999.;
    double tpc_obj_vtx_z = -999.;
    
    double n_pfp{0.0};                // The number of PFP
    int    tpc_obj_mode{0};           // Mode of interaction e.g nue_cc_qe
    int    leading_shower_index{0};   // Index of the leading shower

    // -------------------------------------------------------------------------
    // Initialise the flash variables for this event 
    void SetFlashVariables(std::vector<double> largest_flash_v);
    // -------------------------------------------------------------------------
    // Initialise the TPC object variabeles
    void SetTPCObjVariables(xsecAna::TPCObjectContainer tpc_obj, 
                            double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                            std::vector<double> fv_boundary_v, bool has_pi0);
    // -------------------------------------------------------------------------
    // As above but for DATA 
    void SetTPCObjVariables(xsecAna::TPCObjectContainer tpc_obj, std::string type); 
    // -------------------------------------------------------------------------
    // Helper function to classify the event category
    std::pair<std::string, int> TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool true_in_tpc, bool has_pi0);
    // -------------------------------------------------------------------------
    // Helper function to total up the passed neutrinos and their classifcations
    void TabulateOrigins(std::vector<double> &tabulated_origins, std::string type);
    // -------------------------------------------------------------------------
    // Prints the results of the cuts
    void PrintInfo(int mc_nue_cc_counter, std::vector<double> counter_v, int counter_intime_cosmics,
                              double intime_scale_factor, double data_scale_factor,
                              int counter_dirt, double dirt_scale_factor, std::string cut_name);
    // -------------------------------------------------------------------------
    // Prints total selected in data
    void PrintInfoData(int counter, std::string cut_name);
    // *************************************************************************
    // ------------------- Selection Cuts Functions ----------------------------
    // *************************************************************************
    // Flash is in time and flash PE 
    bool FlashinTime_FlashPE(double flash_time_start, double flash_time_end, double flash_pe_threshold, std::vector<double> &opt_time_v, std::vector<int> &opt_pe_v, std::string type);
    // -------------------------------------------------------------------------
    // Pandora Reco Nue
    bool HasNue(xsecAna::TPCObjectContainer tpc_obj);
    // -------------------------------------------------------------------------
    // In FV
    bool in_fv(double x, double y, double z, std::vector<double> fv_boundary_v);
    // -------------------------------------------------------------------------
    // Flash to vtx distance
    bool opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance);
    bool flashRecoVtxDist(std::vector< double > largest_flash_v, double tolerance, const double tpc_vtx_x, const double tpc_vtx_y, const double tpc_vtx_z);
    // -------------------------------------------------------------------------
    // Vertex to nu distance cut
    bool VtxNuDistance(xsecAna::TPCObjectContainer tpc_obj,int pfp_pdg_type , double tolerance);
    // -------------------------------------------------------------------------
    // Hit Thresholds
    bool HitThreshold(xsecAna::TPCObjectContainer tpc_obj, double threshold, bool useCollection);
    // -------------------------------------------------------------------------
    // Leading Shower Collection Plane Hit Threshold
    bool LeadingHitThreshold(xsecAna::TPCObjectContainer tpc_obj, double threshold);
    // -------------------------------------------------------------------------
    // Open Angle
    bool OpenAngleCut(xsecAna::TPCObjectContainer tpc_obj, double tolerance_open_angle_min, double tolerance_open_angle_max);
    // -------------------------------------------------------------------------
    // dEdx
    bool dEdxCut( xsecAna::TPCObjectContainer tpc_obj, const double tolerance_dedx_min, const double tolerance_dedx_max, std::string type);
    // -------------------------------------------------------------------------
    // Secondary Shower Distance
    bool SecondaryShowersDistCut(xsecAna::TPCObjectContainer tpc_obj, const double dist_tolerance);
    // -------------------------------------------------------------------------
    // Hit per length Ratio
    bool HitLengthRatioCut(const double pfp_hits_length_tolerance, xsecAna::TPCObjectContainer tpc_obj);
    // -------------------------------------------------------------------------
    // Longest Track Leading Shower
    bool LongestTrackLeadingShowerCut(const double ratio_tolerance, xsecAna::TPCObjectContainer tpc_obj);
    // -------------------------------------------------------------------------
    // Contained Tracks Cut
    bool IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v);
    bool ContainedTracksCut(std::vector<double> fv_boundary_v, xsecAna::TPCObjectContainer tpc_obj);
    // -------------------------------------------------------------------------


}; // End Class Selection Cuts

#endif
