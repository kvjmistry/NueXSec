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

    // Set threshold to integrate the flux from [GeV]
    // Remember to make sure you set this number past the bin boundary or it wont work
    // Current flux bins are 0.00 ,0.06, 0.125, 0.25, 0.5... GeV
    double energy_threshold = 0.130; 

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
    
    // Define histograms for the cross section calculation
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> h_cross_sec; // Label -- Universe -- variable -- xsec_type

    // Define histograms for the reweighting by cut
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> h_cut_v; // Label -- Cut -- variable -- Universe

    // enum for histogram vars
    enum TH1D_xsec_hist_vars {
        k_xsec_sel,     // Selected event histogram binned in energy
        k_xsec_bkg,     // Bkg event histogram binned in energy
        k_xsec_gen,     // Gen event histogram binned in energy
        k_xsec_sig,     // Sig event histogram binned in energy
        k_xsec_eff,     // Efficiency histogram binned in energy
        k_xsec_ext,     // EXT event histogram binned in energy
        k_xsec_dirt,    // Dirt event histogram binned in energy
        k_xsec_data,    // Data event histogram binned in energy
        k_xsec_mcxsec,  // MC Cross Section
        k_xsec_dataxsec,// Data Cross Section
        k_TH1D_xsec_MAX
    };

    // enum for histogram vars
    enum TH1D_xsec_var_vars {
        k_var_integrated,     // Integrated X-Section
        k_var_reco_el_E,      // Reconstructed electron energy
        k_var_true_el_E,      // True electron energy
        k_TH1D_xsec_var_MAX
    };

    // Names for cross section histograms
    std::vector<std::string> xsec_types = {"sel", "bkg", "gen", "sig", "eff", "ext", "dirt", "data", "mc_xsec", "data_xsec"};

    std::vector<std::string> vars = {"integrated", "reco_el_E", "true_el_E"
                                    //  "true_nu_E", "reco_nu_e"
                                     };
    
    // Use these for when we do the cross-section
    // std::vector<std::string> var_labels = {";;#nu_{e} + #bar{#nu}_{e} CC Cross-Section [10^{-39} cm^{2}]",
    //                                     ";Reco Leading Shower Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{reco}_{e}} CC Cross-Section [10^{-39} cm^{2}/GeV]",
    //                                     ";True Electron Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{true}_{e}} CC Cross-Section [10^{-39} cm^{2}/GeV]"
    //                                     };
    
    // Use these for when we do the flux normalised event rate
    std::vector<std::string> var_labels = {";;#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}]",
                                        ";Reco. Leading Shower Energy [GeV];#nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}/GeV]",
                                        ";True Electron Energy [GeV]; #nu_{e} + #bar{#nu}_{e} CC Flux Norm. Event Rate [cm^{2}/GeV]"
                                        };
    

    std::vector<std::string> reweighter_labels = {
        "CV",    // Dont comment this out
        "Horn_p2kA",
        "Horn_m2kA",
        "Horn1_x_p3mm",
        "Horm1_x_m3mm",
        "Horn1_y_p3mm",
        "Horn1_y_m3mm",
        "Beam_spot_1_1mm",
        "Beam_spot_1_5mm",
        "Horn2_x_p3mm",
        "Horm2_x_m3mm",
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
        "Dirtup",
        "Dirtdn",
        "POTup",
        "POTdn",
        "pi0",
        "weightsGenie",
        "weightsReint",
        "weightsPPFX",
        "MCStats"
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
        {"Horm1_x_m3mm",        4},
        {"Horn1_y_p3mm",        5},
        {"Horn1_y_m3mm",        6},
        {"Beam_spot_1_1mm",     7},
        {"Beam_spot_1_5mm",     8},
        {"Horn2_x_p3mm",        9},
        {"Horm2_x_m3mm",        10},
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
        k_Horm1_x_m3mm,
        k_Horn1_y_p3mm,
        k_Horn1_y_m3mm,
        k_Beam_spot_1_1mm,
        k_Beam_spot_1_5mm,
        k_Horn2_x_p3mm,
        k_Horm2_x_m3mm,
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
    bool ApplyCuts(int type, SliceContainer &SC, SelectionCuts _scuts);
    // -------------------------------------------------------------------------
    // Fill histograms for each variable for all universes
    void FillCutHists(int type, SliceContainer &SC, std::pair<std::string, int> classification, int cut_index);
    // -------------------------------------------------------------------------
    // Set the weight for universe i depending on the variation 
    void SetUniverseWeight(std::string label, double &weight_uni, double &weight_dirt, double &weight_ext,  double _weightSplineTimesTune,
                           std::string _classification, double cv_weight, int uni, int _nu_pdg, double _true_energy, double _numi_ang, int _npi0, double _pi0_e);
    // -------------------------------------------------------------------------
    // Function to calculate the cross section
    double CalcCrossSec(double sel, double gen, double sig, double bkg, double flux, double ext, double dirt, double targ);
    // -------------------------------------------------------------------------
    // Function to calculate the cross section using binned histograms
    void CalcCrossSecHist(TH1D* h_sel, TH1D* h_eff, TH1D* h_bkg, double mc_scale_factor, double flux, double ext_scale_factor, TH1D* h_ext, double dirt_scale_factor, TH1D* h_dirt, TH1D* h_xsec, double targ, std::string mcdata);
    // -------------------------------------------------------------------------
    // Function to get the integrated flux OR a weight
    double GetIntegratedFlux(int uni, std::string value, std::string label, std::string variation, int _nu_pdg, double _true_energy, double _numi_ang);
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





    private:

    // Here we create the trees 



}; // End Class CrossSectionHelper

#endif