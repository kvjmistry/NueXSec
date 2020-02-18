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


    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type);
    // -------------------------------------------------------------------------
    // Initialise histograms
    void InitHistograms();
    // -------------------------------------------------------------------------
    // Function to make the directory hirarchy
    void MakeDirectory();
    // -------------------------------------------------------------------------
    // Function to fill the reco variables
    void FillReco(int classification_index, int cut_index, SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Function to write the histograms to a file
    void WriteReco(int type);
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
   
    // RecoPlots
    std::vector<std::vector<TH1D*>> h_reco_vtx_x; // Reco Vertex X
    std::vector<std::vector<TH1D*>> h_reco_vtx_y; // Reco Vertex Y
    std::vector<std::vector<TH1D*>> h_reco_vtx_z; // Reco Vertex Z
    std::vector<TH2D*> h_reco_vtx_zy;             // Reco Vertex ZY Plane

    // Leading Shower Momentum
    std::vector<std::vector<TH1D*>> h_reco_leading_mom;

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