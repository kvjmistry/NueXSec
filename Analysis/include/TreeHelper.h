#ifndef TREEHELPER_h
#define TREEHELPER_h

#include "SliceContainer.h"

// Class for filling and saving ttree variables to a file
// This will be useful for producing a flat tree for weighting events
class TreeHelper{

    public:
    // Default constructor
    TreeHelper(){};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type, const char *run_period, const char * file_out );
    // -------------------------------------------------------------------------
    // Function to fill the tree vars
    void FillVars(SliceContainer &SC, std::pair<std::string, int> _classification, bool _gen, double _weight, double _reco_energy);
    // -------------------------------------------------------------------------
    // Fill the variables in the dedx tree
    void Fill_dedxVars(SliceContainer &SC, std::pair<std::string, int> _classification, std::string _cut, double _weight);
    // -------------------------------------------------------------------------
    // Writes the tree to file
    void WriteTree();
    // -------------------------------------------------------------------------
    // Fill the counter tree
    void Fill_counters(std::vector<double> counter_v, bool bool_use_mc, bool bool_use_ext, bool bool_use_data, bool bool_use_dirt);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    
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

    double true_energy{0.0}, reco_energy{0.0};
    int n_showers{0}, n_tracks{0};
    float shr_phi{0.0};
    float shr_energy_tot_cali{0.0};
    float shrmoliereavg{0.0};
    float shr_hits_max{0.0};
    float elec_e{0.0};

    // Weights
    std::vector<unsigned short> weightsGenie;
    std::vector<unsigned short> weightsReint;
    std::vector<unsigned short> weightsPPFX ;
    double knobRPAup;
    double knobCCMECup;
    double knobAxFFCCQEup;
    double knobVecFFCCQEup;
    double knobDecayAngMECup;
    double knobThetaDelta2Npiup;
    double knobThetaDelta2NRadup;
    double knobRPA_CCQE_Reducedup;
    double knobNormCCCOHup;
    double knobNormNCCOHup;
    double knobRPAdn;
    double knobCCMECdn;
    double knobAxFFCCQEdn;
    double knobVecFFCCQEdn;
    double knobDecayAngMECdn;
    double knobThetaDelta2Npidn;
    double knobThetaDelta2NRaddn;
    double knobRPA_CCQE_Reduceddn;
    double knobNormCCCOHdn;
    double knobNormNCCOHdn;



    TTree * tree;       // Main tree with the selected events
    TTree * dedx_tree;  // Tree for optimising the dedx cut
    TTree * counter_tree; // Tree for storing the selection results


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

    // Counters
    double count_nue_cc_qe{0.0};
    double count_nue_cc_res{0.0};
    double count_nue_cc_dis{0.0};
    double count_nue_cc_coh{0.0};
    double count_nue_cc_mec{0.0};
    
    double count_nuebar_cc_qe{0.0};
    double count_nuebar_cc_res{0.0};
    double count_nuebar_cc_dis{0.0};
    double count_nuebar_cc_coh{0.0};
    double count_nuebar_cc_mec{0.0};
    
    double count_nue_cc_infv{0.0};
    double count_nuebar_cc_infv{0.0};
    double count_nue_cc_incryo{0.0};
    double count_nuebar_cc_incryo{0.0};
    
    double count_numu_cc_qe{0.0};
    double count_numu_cc_res{0.0};
    double count_numu_cc_dis{0.0};
    double count_numu_cc_coh{0.0};
    double count_numu_cc_mec{0.0};
    
    double count_numubar_cc_qe{0.0};
    double count_numubar_cc_res{0.0};
    double count_numubar_cc_dis{0.0};
    double count_numubar_cc_coh{0.0};
    double count_numubar_cc_mec{0.0};
    
    double count_numu_cc_infv{0.0};
    double count_numubar_cc_infv{0.0};
    double count_numu_cc_incryo{0.0};
    double count_numubar_cc_incryo{0.0};
    
    double count_nue_cc{0.0};
    double count_nuebar_cc{0.0};
    double count_nu_out_fv{0.0};
    double count_cosmic{0.0};
    double count_numu_cc{0.0};
    double count_numu_cc_pi0{0.0};
    double count_nc{0.0};
    double count_nc_pi0{0.0};
    double count_unmatched{0.0};
    double count_total_mc{0.0};
    double count_data{0.0};
    double count_ext{0.0};
    double count_dirt{0.0};

    

    private:

    // Here we create the trees 



}; // End Class TreeHelper

#endif