#ifndef CROSSSECTIONHELPER_H
#define CROSSSECTIONHELPER_H

#include "Utility.h"
#include "SliceContainer.h"
#include "SelectionCuts.h"

// Class for calculating the cross section and doing systematics
class CrossSectionHelper{

    public:
    // Default constructor
    CrossSectionHelper(){};
    
    // The input and output files
    TFile* f_nuexsec, *fnuexsec_out, *f_flux;

    // Name of the merged flux file
    std::string flux_file_name = "";

    // Class instances
    Utility _util;

    // Variables
    int run{0}, subrun{0}, event{0};
    std::string *classification = NULL; // The classification of the event
    
    // Is the event signal?
    bool gen{false};

    // Did the event pass the selection? (some signal events in the eff denominator wont)
    bool passed_selection{false};

    double weight{0.0};        // This is not going to be integer if we already weight the CV

    double true_energy{0.0}, reco_energy{0.0};
    float shr_energy_cali{0.0};
    float elec_e{0.0};
    float ppfx_cv{1.0};
    float weightSplineTimesTune{1.0};
    float numi_ang{0.0};
    int nu_pdg{0};
    int npi0{0};
    double pi0_e{0.0};
    double effective_angle{0.0};     // The angle between the vector from the target to nu vtx compared to the reconstructed shower direction.
    double cos_effective_angle{0.0}; // The cosine of the angle between the vector from the target to nu vtx compared to the reconstructed shower direction.
    double true_effective_angle{0.0};     // True angle between electron and neutrino vectors
    double cos_true_effective_angle{0.0}; // True angle between electron and neutrino vectors
    int ccnc{0};

    // Weights
    std::vector<unsigned short> *weightsGenie = NULL;
    std::vector<unsigned short> *weightsReint = NULL;
    std::vector<unsigned short> *weightsPPFX = NULL ;
    double knobRPAup{0.0};
    double knobCCMECup{0.0};
    double knobAxFFCCQEup{0.0};
    double knobVecFFCCQEup{0.0};
    double knobDecayAngMECup{0.0};
    double knobThetaDelta2Npiup{0.0};
    double knobThetaDelta2NRadup{0.0};
    double knobRPA_CCQE_Reducedup{0.0};
    double knobNormCCCOHup{0.0};
    double knobNormNCCOHup{0.0};
    double knobxsr_scc_Fv3up{0.0};
    double knobxsr_scc_Fa3up{0.0};
    double knobRPAdn{0.0};
    double knobCCMECdn{0.0};
    double knobAxFFCCQEdn{0.0};
    double knobVecFFCCQEdn{0.0};
    double knobDecayAngMECdn{0.0};
    double knobThetaDelta2Npidn{0.0};
    double knobThetaDelta2NRaddn{0.0};
    double knobRPA_CCQE_Reduceddn{0.0};
    double knobNormCCCOHdn{0.0};
    double knobNormNCCOHdn{0.0};
    double knobxsr_scc_Fv3dn{0.0};
    double knobxsr_scc_Fa3dn{0.0};

    std::vector<double> vec_universes;

    TTree * tree;

    std::string run_period;

    // Cross section Variables -- So far copied from coltons analysis, so these numbers need updating
    double lar_density_mc   = 1.3836;   // Density of Argon in the simulation g/cm3
    double lar_density_data = 1.3836;   // Density of Argon in the simulation g/cm3
    double volume           = -1; // Fiducial volume cm3
    double NA               = 6.022140857e23; // Advogadro's number molecule/mol
    double N_nuc            = 40.0;     // Number of argon nuclons
    double m_mol            = 39.95;    // Molar Mass of Argon g/mol

    double N_target_MC{0.0};   // Set in code, total number of targets in MC
    double N_target_Data{0.0}; // Set in code, total number of targets in Data

    // Fluxes need to be defined by reading in the flux file and integrating
    double integrated_flux{0.0};
    double flux_scale_factor{1.0e-4}; // unit conversion of flux from m2 to cm2
    double mc_flux_scale_factor{1.0};
    double data_flux_scale_factor{1.0};

    int uni_reint{1000}, uni_genie{600}, uni_ppfx{600}, uni_mcstats{1000}; // For resizing data, ext and dirt in multisims

    // Random number generator for generating MC stats uncertainty
    TRandom3 *rand = new TRandom3(); 

    // Define the bins for each variable -- See InitialiseHistograms function to see the actual binning used
    std::vector<std::vector<double>> bins;
    std::vector<std::vector<double>> bins_fine;

    // Define histograms for the cross section calculation
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> h_cross_sec; // Label -- Universe -- variable -- xsec_type

    // define vector of smearing matrices
    std::vector<std::vector<std::vector<TH2D*>>> h_smear; // Label -- Universe -- variable
    std::vector<std::vector<std::vector<TH2D*>>> h_smear_fine; // Label -- Universe -- variable

    // Define histograms for the reweighting by cut
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> h_cut_v; // Label -- Cut -- variable -- Universe

    TH2D* h_2D_CV_binning; // Histogram to define the bin indexes for a response matrix in energy and angle

    TH1D *h_detvar_cv; // the detector variation cv for us to apply the response matrix to. 

    // enum for histogram types
    enum TH1D_xsec_hist_vars {
        k_xsec_sel,          // Selected event histogram binned in energy
        k_xsec_bkg,          // Bkg event histogram binned in energy
        k_xsec_gen,          // Gen event histogram binned in energy
        k_xsec_gen_smear,    // Gen event histogram binned in energy with smeared truth
        k_xsec_gen_shape,    // Generated events for shape uncertainties
        k_xsec_sig,          // Sig event histogram binned in energy
        k_xsec_eff,          // Efficiency histogram binned in energy
        k_xsec_ext,          // EXT event histogram binned in energy
        k_xsec_dirt,         // Dirt event histogram binned in energy
        k_xsec_data,         // Data event histogram binned in energy
        k_xsec_mcxsec,       // MC Cross Section
        k_xsec_mcxsec_smear, // MC Cross Section smeared truth
        k_xsec_mcxsec_shape, // MC Cross Section shape uncertainties
        k_xsec_dataxsec,     // Data Cross Section
        k_TH1D_xsec_MAX
    };

    // enum for histogram vars
    enum TH1D_xsec_var_vars {
        k_var_integrated,     // Total X-Section
        k_var_recoX,          // X-Sec as a function of a Reconstructed variable
        k_var_trueX,          // X-Sec as a function of a True variable
        k_TH1D_xsec_var_MAX
    };

    // Names for cross section histograms
    std::vector<std::string> xsec_types = {"sel", "bkg", "gen", "gen_smear", "gen_shape", "sig", "eff", "ext", "dirt", "data", "mc_xsec", "mc_xsec_smear", "mc_xsec_shape", "data_xsec"};

    std::vector<std::string> vars = {"integrated", "recoX", "trueX"};
    
    // Use these for when we do the flux normalised event rate
    // std::vector<std::string> var_labels = {";;#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}]",
    //                                     ";Reco. Leading Shower Energy [GeV];#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}/GeV]",
    //                                     ";True Electron Energy [GeV]; #nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}/GeV]"
    //                                     };
    
    std::vector<std::string> var_labels_xsec = {};

    std::vector<std::string> var_labels_events = {};

    std::vector<std::string> var_labels_eff = {};

    std::string smear_hist_name = ";True e#lower[-0.5]{-} + e^{+} Energy [GeV];Leading Shower Energy [GeV]";

    std::vector<double> hist_bins;

    std::vector<double> fine_bins;

    double xsec_scale = 0.0; // not needed in this code

    std::vector<std::string> reweighter_labels = {
        "CV",    // Dont comment this out
        "Horn_p2kA",
        "Horn_m2kA",
        "Horn1_x_p3mm",
        "Horn1_x_m3mm",
        "Horn1_y_p3mm",
        "Horn1_y_m3mm",
        "Beam_spot_1_1mm",
        "Beam_spot_1_5mm",
        "Horn2_x_p3mm",
        "Horn2_x_m3mm",
        "Horn2_y_p3mm",
        "Horn2_y_m3mm",
        "Horns_0mm_water",
        "Horns_2mm_water",
        "Beam_shift_x_p1mm",
        "Beam_shift_x_m1mm",
        "Beam_shift_y_p1mm",
        "Beam_shift_y_m1mm",
        "Target_z_p7mm",
        "Target_z_m7mm",
        "Horn1_refined_descr",
        "Decay_pipe_Bfield",
        "Old_Horn_Geometry",
        "RPAup",
        "CCMECup",
        "AxFFCCQEup",
        "VecFFCCQEup",
        "DecayAngMECup",
        "ThetaDelta2Npiup",
        "ThetaDelta2NRadup",
        "RPA_CCQE_Reducedup",
        "NormCCCOHup",
        "NormNCCOHup",
        "RPAdn",
        "CCMECdn",
        "AxFFCCQEdn",
        "VecFFCCQEdn",
        "DecayAngMECdn",
        "ThetaDelta2Npidn",
        "ThetaDelta2NRaddn",
        "RPA_CCQE_Reduceddn",
        "NormCCCOHdn",
        "NormNCCOHdn",
        "xsr_scc_Fv3up",
        "xsr_scc_Fa3up",
        "xsr_scc_Fv3dn",
        "xsr_scc_Fa3dn",
        "Dirtup",
        "Dirtdn",
        "POTup",
        "POTdn",
        "pi0",
        "weightsGenie",
        "weightsReint",
        "weightsPPFX"
        // "MCStats"
    };

    std::vector<std::vector<TH2D*>> beamline_hists;
    std::vector<double> beamline_flux;

    enum nu_flav {
        k_numu,
        k_numubar,
        k_nue,
        k_nuebar,
        k_flav_max
    };

    // vector fo going to and from each beamline pair
    std::vector< std::pair<std::string,int>> beamline_map= {
        {"Horn_p2kA",           1}, 
        {"Horn_m2kA",           2}, 
        {"Horn1_x_p3mm",        3},
        {"Horn1_x_m3mm",        4},
        {"Horn1_y_p3mm",        5},
        {"Horn1_y_m3mm",        6},
        {"Beam_spot_1_1mm",     7},
        {"Beam_spot_1_5mm",     8},
        {"Horn2_x_p3mm",        9},
        {"Horn2_x_m3mm",        10},
        {"Horn2_y_p3mm",        11},
        {"Horn2_y_m3mm",        12},
        {"Horns_0mm_water",     13},
        {"Horns_2mm_water",     14},
        {"Beam_shift_x_p1mm",   15},
        {"Beam_shift_x_m1mm",   16},
        {"Beam_shift_y_p1mm",   17},
        {"Beam_shift_y_m1mm",   18},
        {"Target_z_p7mm",       19},
        {"Target_z_m7mm",       20},
        {"Horn1_refined_descr", 21},
        {"Decay_pipe_Bfield",   22},
        {"Old_Horn_Geometry",   23}
    };

    enum beamline_enums {
        k_Horn_p2kA,
        k_Horn_m2kA,
        k_Horn1_x_p3mm,
        k_Horn1_x_m3mm,
        k_Horn1_y_p3mm,
        k_Horn1_y_m3mm,
        k_Beam_spot_1_1mm,
        k_Beam_spot_1_5mm,
        k_Horn2_x_p3mm,
        k_Horn2_x_m3mm,
        k_Horn2_y_p3mm,
        k_Horn2_y_m3mm,
        k_Horns_0mm_water,
        k_Horns_2mm_water,
        k_Beam_shift_x_p1mm,
        k_Beam_shift_x_m1mm,
        k_Beam_shift_y_p1mm,
        k_Beam_shift_y_m1mm,
        k_Target_z_p7mm,
        k_Target_z_m7mm,
        k_Horn1_refined_descr,
        k_Decay_pipe_Bfield,
        k_Old_Horn_Geometry,
        k_beamline_max
    };


    // Filelists
    // Create a file for the selected events and their properties
    std::ofstream evt_dist_sig, evt_dist_gen, evt_dist_bkg, evt_dist_data;

    // Vector for storing event weights
    std::vector<float> ev_weight;

    // To initialise file names for event list readout
    bool filled_sig = false, filled_bkg = false, filled_gen = false;

    // For Filling TTree in Andy's analysis format
    TFile *f_out;
    TTree *event_tree;
    TTree *meta_tree;

    float potData = 2.0;
    float potMC   = 23.2136;
    float n_targ  = 4.31247;
    float nUniverses = 600;
    bool isData, isSignal, isSelected; 
    float xTrue, xReco;
    double integrated_flux_tree{0.0};

    double recoX{0.0};
    double trueX{0.0};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(Utility _utility);
    // -------------------------------------------------------------------------
    // Function to loop over events and calculate the cross section
    void LoopEvents(); 
    // -------------------------------------------------------------------------
    // Loop over events from initial tree and reweight them to get the sys uncertainty
    // by cut type for a number of variables
    void LoopEventsbyCut();
    // -------------------------------------------------------------------------
    // Apply the selection cuts -- so we can reweight events at differnt points in the selection
    bool ApplyCuts(int type, SliceContainer &SC, SelectionCuts _scuts, int treeNum);
    // -------------------------------------------------------------------------
    // Fill histograms for each variable for all universes
    void FillCutHists(int type, SliceContainer &SC, std::pair<std::string, int> classification, int cut_index);
    // -------------------------------------------------------------------------
    // Set the weight for universe i depending on the variation 
    void SetUniverseWeight(std::string label, double &weight_uni, double &weight_dirt, double _weightSplineTimesTune,
                           std::string _classification, double cv_weight, int uni, int _nu_pdg, double _true_energy, double _numi_ang, int _npi0, double _pi0_e, int _ccnc);
    // -------------------------------------------------------------------------
    // Function to calculate the cross section
    double CalcCrossSec(double sel, double gen, double sig, double bkg, double flux, double ext, double dirt, double targ);
    // -------------------------------------------------------------------------
    // Function to calculate the cross section using binned histograms
    void CalcCrossSecHist(TH1D* h_sel, TH1D* h_eff, TH1D* h_bkg, double mc_scale_factor, double flux, double ext_scale_factor, TH1D* h_ext, double dirt_scale_factor, TH1D* h_dirt, TH1D* h_xsec, TH1D* h_sig, double targ, std::string mcdata, int _var);
    // -------------------------------------------------------------------------
    // Function to get the integrated flux for CV
    double GetIntegratedFluxCV();
    // -------------------------------------------------------------------------
    // Get the integrated flux from the FLUGG flux file
    double GetIntegratedFluxFLUGG();
    // -------------------------------------------------------------------------
    // Function to get the integrated flux for HP Universe
    double GetIntegratedFluxHP(int uni, std::string label);
    // -------------------------------------------------------------------------
    // Function to get the weight for a beamline variation
    double GetBeamlineWeight(std::string variation, int _nu_pdg, double _true_energy, double _numi_ang);
    // -------------------------------------------------------------------------
    // Function to get the weight for a HP variation using the flux file
    double GetHPWeight(int uni, std::string label, int _nu_pdg, double _true_energy, double _numi_ang);
    // -------------------------------------------------------------------------
    // Function to get the POT from the flux file
    double GetPOT(TFile* f);
    // -------------------------------------------------------------------------
    int GetBinIndex(double reco_energy);
    // -------------------------------------------------------------------------
    void WriteHists();
    // -------------------------------------------------------------------------
    // Initialise the ttree reading from the input file
    void InitTree();
    // -------------------------------------------------------------------------
    // Function will set the reweight vector to the corresponding label per event
    void SwitchReweighterLabel(std::string label);
    // -------------------------------------------------------------------------
    // Override the above functio, but use the slice container class that has been initialised
    void SwitchReweighterLabel(std::string label, SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Initialise the histograms for this class
    void InitialiseHistograms(std::string run_mode);
    // -------------------------------------------------------------------------
    // Helper function to fill the histograms
    void FillHists(int label, int uni, int xsec_type, double weight_uni, float shr_energy_cali, float elec_e);
    // -------------------------------------------------------------------------
    // Function to get the beamline files and open them
    void GetBeamlineHists();
    // -------------------------------------------------------------------------
    // Checks whether the variation is a beamline
    // We use this to override the weights and integrated flux so we can get them
    // from the beamline file
    bool CheckBeamline(std::string variation);
    // -------------------------------------------------------------------------
    // Get the index of the beamline variation so we can get the correct histogram to weight from
    int GetBeamlineIndex(std::string variation);
    // -------------------------------------------------------------------------
    // Generates a Poisson random number
    double PoissonRandomNumber(int uni);
    // -------------------------------------------------------------------------
    // Concatenates the run subrun event number to a single number
    int ConcatRunSubRunEvent(int run, int subrun, int event);
    // -------------------------------------------------------------------------
    // Normalise the smearing matrix and then use it to get a smeared efficiency
    void Smear(TH1D* h_sig, TH1D* h_gen, TH2D* h_smear, TH1D* h_eff);
    // -------------------------------------------------------------------------
    // Create a response matrix and smear the generated events in truth to sig events in reco
    void ApplyResponseMatrix(TH1D* h_gen, TH1D* h_gen_smear, TH1D* h_gen_CV, TH2D* h_smear, bool norm);
    // -------------------------------------------------------------------------
    // Save the event to file
    void SaveEvent(std::string _classification, bool _passed_selection, std::vector<float> ev_weight, double reco_E, double true_E);
    // -------------------------------------------------------------------------
    // Return the bin index from a 2D histogram
    int GetBinIndex();
    // -------------------------------------------------------------------------
    // Load in the detector variation CV so we can smear it
    void LoadDetvarCVHist();
    // -------------------------------------------------------------------------
    // Unfold the data in an unregularized way
    void UnregularizedUnfold(TH1D *h_data_xsec_reco, TH1D* h_data_xsec_true, TH2D* h_response);
    // -------------------------------------------------------------------------
    // Make and save generator histograms to file
    void SaveGenXSec();
    // -------------------------------------------------------------------------
    // Compare the nuwro and genie pi0 cross sections
    void CheckPi0CrossSection();






    private:

    // Here we create the trees 



}; // End Class CrossSectionHelper

#endif