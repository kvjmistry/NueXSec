#ifndef SLICECONTAINER_H
#define SLICECONTAINER_H

#include "utility.h"

/* 
Class to hold information for the eventslice in the pandora ntuples for ease of use

*/

// Slice Container Class
class SliceContainer {
public:

    // -------------------------------------------------------------------------
    // Initialise the class
    void Initialise(TTree *tree, TTree *mc_truth_tree, int type, TFile *f_flux_weights, const char * _run_period);
    // -------------------------------------------------------------------------
    // Function to classify the slice
    std::pair<std::string, int>  SliceClassifier(int type);
    // -------------------------------------------------------------------------
    // Function to return the genie interaction mode, e.g. ccqe, ccmec etc.
    std::string SliceInteractionType(int type);
    // -------------------------------------------------------------------------
    // Function to classify the event by particle type of the leading shower
    std::pair<std::string, int> ParticleClassifier(int type);
    // -------------------------------------------------------------------------
    // Reset the particle container objects
    void Reset();
    // -------------------------------------------------------------------------
    // Get the Leading shower vertex
    void GetLeadingShowerIndex();
    // -------------------------------------------------------------------------
    // Check if true neutrino vertex is in the FV
    bool InFV();
    // -------------------------------------------------------------------------
    // Set the TPC object particle tree
    void SetTPCObj(xsecAna::TPCObjectContainer &tpc_obj, int type);
    // -------------------------------------------------------------------------


    utility _util;


    // TPC Object Container
    std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;

    int n_pfp;

    int longest_track_index;
    int leading_shower_index;
    int subleading_shower_index;
    int n_pfp_tracks;
    int n_pfp_showers;
    int tpc_obj_mode;

    double tpco_vtx_x;
    double tpco_vtx_y;
    double tpco_vtx_z;
    double tpco_hits;

    bool infv; // True nu vtx in the FV

    int run;
    int subrun;
    int event;

    std::string tpc_obj_origin;

    // Particle Properties
    std::vector<std::string> mc_origin;
    std::vector<int>  pfp_pdg;
    std::vector<int>  num_pfp_hits;
    std::vector<int>  mc_parent_pdg;
    std::vector<int>  n_pfp_hits_w;
  
    std::vector<double> pfp_length;
    std::vector<double> pfp_open_angle;
    std::vector<double> leading_dedx;
  
    std::vector<double> pfp_vtx_x;
    std::vector<double> pfp_vtx_y;
    std::vector<double> pfp_vtx_z;
    std::vector<double> pfp_dir_x;
    std::vector<double> pfp_dir_y;
    std::vector<double> pfp_dir_z;


    std::vector<double> mc_pdg;
    std::vector<double> mc_Theta;
    std::vector<double> mc_Phi;
    std::vector<double> mc_Energy;
    
    std::vector<double> pfp_end_x;
    std::vector<double> pfp_end_y;
    std::vector<double> pfp_end_z;

    std::vector<double> dist_pfp_nu_vtx;

    std::vector<int> CCNC;

    // MC Truth Counter Tree
    int mc_nue_cc_counter = 0;
    int mc_nue_nc_counter = 0;
    int mc_numu_cc_counter = 0;
    int mc_numu_nc_counter = 0;
    int mc_nue_cc_counter_bar = 0;
    int mc_numu_cc_counter_bar = 0;
    int mc_nue_nc_counter_bar = 0;
    int mc_numu_nc_counter_bar = 0;
    double mc_nu_energy = 0;
    double mc_nu_momentum = 0;
    int mc_nu_id;

    int mc_nu_num_particles = 0;
    int mc_nu_num_charged_particles = 0;
    
    double mc_nu_vtx_x = -999;
    double mc_nu_vtx_y = -999;
    double mc_nu_vtx_z = -999;
    double mc_nu_dir_x = -999;
    double mc_nu_dir_y = -999;
    double mc_nu_dir_z = -999;
    
    double mc_ele_dir_x = -999;
    double mc_ele_dir_y = -999;
    double mc_ele_dir_z = -999;
    
    double mc_ele_energy = 0;
    double mc_ele_momentum = 0;
    bool has_pi0 = false;
    double mc_nu_time = -1;


};

#endif
