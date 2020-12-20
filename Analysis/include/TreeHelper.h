#ifndef TREEHELPER_H
#define TREEHELPER_H

#include "SliceContainer.h"

// Class for filling and saving ttree variables to a file
// This will be useful for producing a flat tree for weighting events
class TreeHelper{

    public:
    // Default constructor
    TreeHelper(){};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type, const char* file_out, Utility _utility);
    // -------------------------------------------------------------------------
    // Function to fill the tree vars
    void FillVars(SliceContainer &SC, bool passed_selection);
    // -------------------------------------------------------------------------
    // Fill the variables in the dedx tree
    void Fill_dedxVars(SliceContainer &SC, std::pair<std::string, int> _classification, std::string _cut, double _weight);
    // -------------------------------------------------------------------------
    // Writes the tree to file
    void WriteTree(int type);
    // -------------------------------------------------------------------------
    // Fill the counter tree
    void Fill_counters(std::vector<double> counter_v, bool bool_use_mc, bool bool_use_ext, bool bool_use_data, bool bool_use_dirt);
    // -------------------------------------------------------------------------
    // Set the branch address of some ttrees
    void SetBranches(TTree * tree);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    
    // Destructor 
    // ~TreeHelper(); 

    // The output file
    TFile* f_nuexsec;

    // Class instances
    Utility _util;

    int _type{1};

    // Tree variables
    int run{0}, subrun{0}, event{0};
    std::string classification; // The classification of the event
    
    // Is the event a true signal event in the FV that was not selected?
    // We still need these for the efficiency
    bool gen{false};   
    bool passed_selection{false};        
    
    double weight{0.0};        // This is not going to be integer if we already weight the CV

    double true_energy{0.0}, reco_energy{0.0};
    int n_showers{0}, n_tracks{0};
    float shr_phi{0.0};
    float shr_energy_cali{0.0};
    float shrmoliereavg{0.0};
    float shr_hits_max{0.0};
    float elec_e{0.0};
    float ppfx_cv{1.0};
    float weightSplineTimesTune{1.0};
    float numi_ang{0};
    int nu_pdg{0};
    int shr_bkt_pdg{0};
    float shr_bkt_purity{0.0};
    float shr_bkt_completeness{0.0};
    float shr_bkt_E{0.0}; // energy of truth matched particle to the leading shower
    int npi0{0};
    double pi0_e{0.0};
    int interaction{0};
    std::vector<float> all_shr_hits;
    std::vector<float> all_shr_energies;


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



    TTree * tree;             // Main tree with the selected events
    TTree * nue_tree;         // Main tree with the selected events for intrinsic nue sample
    TTree * dedx_tree;        // Tree for optimising the dedx cut
    TTree * counter_tree;     // Tree for storing the selection results
    TTree * nue_counter_tree; // Tree for storing the nue selection results


    // vars for dedx
    unsigned int shr_hits_u_tot{0},         shr_hits_v_tot{0},         shr_hits_y_tot{0};
    float shr_tkfit_dedx_Y{0.0},       shr_tkfit_dedx_V{0.0},       shr_tkfit_dedx_U{0.0};
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

    double count_pi0_nue_cc_nopi0{0.0};
    double count_pi0_nue_cc_pi0{0.0};
    double count_pi0_nuebar_cc_nopi0{0.0};
    double count_pi0_nuebar_cc_pi0{0.0};
    double count_pi0_numu_cc_nopi0{0.0};
    double count_pi0_numu_cc_pi0{0.0};
    double count_pi0_nc_nopi0{0.0};
    double count_pi0_nc_pi0{0.0};
    
    double count_nue_cc{0.0};
    double count_nuebar_cc{0.0};
    double count_nu_out_fv{0.0};
    double count_cosmic{0.0};
    double count_numu_cc{0.0};
    double count_numu_cc_pi0{0.0};
    double count_nc{0.0};
    double count_nc_pi0{0.0};
    double count_unmatched{0.0};
    double count_unmatched_nue{0.0};
    double count_cosmic_nue{0.0};
    double count_unmatched_nuebar{0.0};
    double count_cosmic_nuebar{0.0};
    double count_total_mc{0.0};
    double count_data{0.0};
    double count_ext{0.0};
    double count_dirt{0.0};

    

    private:

    // Here we create the trees 



}; // End Class TreeHelper

#endif