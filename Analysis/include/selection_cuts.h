#ifndef SELECTION_CUTS_h
#define SELECTION_CUTS_h

#include "utility.h"

// Class for applying the selection cuts. Cut classes will inherit
// from this one and override function
class selection_cuts{

    private:

    // Bools to switch file type 
    bool mc, data, ext, dirt, variation; 

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
    
    double n_pfp;                // The number of PFP
    int    tpc_obj_mode;         // Mode of interaction e.g nue_cc_qe
    int    leading_shower_index; // Index of the leading shower

    // -------------------------------------------------------------------------
    // Initialise the flash variables for this event 
    void SetFlashVariables(std::vector<double> largest_flash_v);
    // -------------------------------------------------------------------------
    // Initialise the TPC object variabeles
    void SetTPCObjVariables(xsecAna::TPCObjectContainer tpc_obj, 
                            double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                            std::vector<double> fv_boundary_v, bool has_pi0);
    // -------------------------------------------------------------------------
    // Helper function to classify the event category
    std::pair<std::string, int> TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool true_in_tpc, bool has_pi0);
    // -------------------------------------------------------------------------

    // *************************************************************************
    // ------------------- Selection Cuts Functions ----------------------------
    // *************************************************************************
    // Flash is in time and flash PE 
    void FlashinTime_FlashPE(TFile* f, double flash_start_time, double flash_end_time, std::vector<bool> &flash_cuts_pass_vec, TString mode );
    // -------------------------------------------------------------------------
    // Pandora Reco Nue
    bool HasNue(xsecAna::TPCObjectContainer tpc_obj, const int n_pfp );
    // -------------------------------------------------------------------------
    // In FV
    bool in_fv(double x, double y, double z, std::vector<double> fv_boundary_v);
    // -------------------------------------------------------------------------
    // Flash to vtx distance
    bool flashRecoVtxDist(std::vector< double > largest_flash_v, double tolerance, const double tpc_vtx_x, const double tpc_vtx_y, const double tpc_vtx_z);
    // -------------------------------------------------------------------------
    // Vertex to nu distance cut
    bool VtxNuDistance(xsecAna::TPCObjectContainer tpc_obj,int pfp_pdg_type , double tolerance);
    // -------------------------------------------------------------------------
    // Hit Thresholds
    bool HitThreshold(xsecAna::TPCObjectContainer tpc_obj, double threshold, bool useCollection);
    // -------------------------------------------------------------------------
    // Open Angle
    bool OpenAngleCut(xsecAna::TPCObjectContainer tpc_obj, const std::vector<double> tolerance_open_angle);
    // -------------------------------------------------------------------------
    // dEdx
    bool dEdxCut( xsecAna::TPCObjectContainer tpc_obj, const double tolerance_dedx_min, const double tolerance_dedx_max);
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
    bool ContainedTracksCut(std::vector<double> fv_boundary_v, xsecAna::TPCObjectContainer tpc_obj);
    // -------------------------------------------------------------------------


}; // End Class Selection Cuts

#endif
