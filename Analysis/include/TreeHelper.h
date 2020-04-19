#ifndef TREEHELPER_h
#define TREEHELPER_h

#include "SliceContainer.h"

// Class for filling and saving ttree variables to a file
// This will be useful for producing a flat tree for weighting events
class TreeHelper{

    public:
    // Default constructor
    TreeHelper(){};
    
    // Destructor 
    // ~TreeHelper(); 

    // The output file
    TFile* f_nuexsec;

    // Class instances
    utility _util;

    int _type{1};

    bool weight_tune = true; // Apply genie tune weight
    bool weight_ppfx = true; // Apply ppfx cv weight

    // Tree variables
    int run{0}, subrun{0}, event{0};
    std::string classifcation; // The classification of the event
    
    // Is the event a true signal event in the FV that was not selected?
    // We still need these for the efficiency
    bool gen{false};           
    
    double weight{0.0};        // This is not going to be integer if we already weight the CV

    double efficiency{0.0}, purity{0.0};

    TTree * tree;     // Main tree with the selected events
    TTree * eff_tree; // Efficiency and Purity tree
    TTree * dedx_tree;  // Tree for optimising the dedx cut


    // vars for dedx
    float shr_dedx_Y_cali{0.0},        shr_dedx_V_cali{0.0},        shr_dedx_U_cali{0.0};
    float shr_tkfit_dedx_Y{0.0},       shr_tkfit_dedx_V{0.0},       shr_tkfit_dedx_U{0.0};
    float shr_tkfit_dedx_Y_alt{0.0},   shr_tkfit_dedx_V_alt{0.0},   shr_tkfit_dedx_U_alt{0.0};
    float shr_tkfit_2cm_dedx_Y{0.0},   shr_tkfit_2cm_dedx_V{0.0},   shr_tkfit_2cm_dedx_U{0.0};
    float shr_tkfit_gap05_dedx_Y{0.0}, shr_tkfit_gap05_dedx_V{0.0}, shr_tkfit_gap05_dedx_U{0.0};
    float shr_tkfit_gap10_dedx_Y{0.0}, shr_tkfit_gap10_dedx_V{0.0}, shr_tkfit_gap10_dedx_U{0.0};
    float shr_distance{0.0};
    float shr_theta{0.0};

    std::string cut;

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type, const char *run_period, std::string file_out, int weight_cfg );
    // -------------------------------------------------------------------------
    // Function to fill the tree vars
    void FillVars(SliceContainer &SC, std::pair<std::string, int> _classification, bool _gen, double _weight);
    // -------------------------------------------------------------------------
    // Fill the variables in the dedx tree
    void Fill_dedxVars(SliceContainer &SC, std::pair<std::string, int> _classification, std::string _cut, double _weight);
    // -------------------------------------------------------------------------
    // Fills the Efficiency and Purity
    void FillEff(double _efficiency, double _purity);
    // -------------------------------------------------------------------------
    // Writes the tree to file
    void WriteTree(int type);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    private:

    // Here we create the trees 



}; // End Class TreeHelper

#endif