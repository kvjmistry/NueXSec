#ifndef CROSSSECTIONHELPER_H
#define CROSSSECTIONHELPER_H

#include "Utility.h"

// Class for calculating the cross section and doing systematics
class CrossSectionHelper{

    public:
    // Default constructor
    CrossSectionHelper(){};
    
    // The input and output files
    TFile* f_nuexsec, *fnuexsec_out;

    // Class instances
    Utility _util;

     // Scale factors (everything is scaled to data)
    double mc_scale_factor     = 1.0;
    double ext_scale_factor     = 1.0;
    double dirt_scale_factor   = 1.0;

    // Variables
    int run{0}, subrun{0}, event{0};
    std::string *classifcation = NULL; // The classification of the event
    
    // Is the event a true signal event in the FV that was not selected?
    // We still need these for the efficiency
    bool gen{false};           
    
    double weight{0.0};        // This is not going to be integer if we already weight the CV

    double true_energy{0.0}, reco_energy{0.0};
    float shr_energy_cali{0.0};
    float elec_e{0.0};
    float ppfx_cv{1.0};
    float weightSplineTimesTune{1.0};
    float numi_ang{0.0};
    int nu_pdg{0};

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
    double lar_density_mc   = 1.3954;   // Density of Argon in the simulation g/cm3
    double lar_density_data = 1.3836;   // Density of Argon in the simulation g/cm3
    double volume           = 4.1622e7; // Fiducial volume cm3
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

    int uni_reint{1000}, uni_genie{600}, uni_ppfx{600}; // For resizing data, ext and dirt in multisims

    // Define the bins for each variable -- See InitialiseHistograms function to see the actual binning used
    std::vector<std::vector<double>> bins;
    
    // Define histograms for the cross section calculation
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> h_cross_sec; // Label -- Universe -- variable -- xsec_type

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
        // k_var_true_nu_E,      // True neutrino energy
        // k_var_reco_nu_E,      // Reconstructed neutrino energy
        k_TH1D_xsec_var_MAX
    };

    // Names for cross section histograms
    std::vector<std::string> xsec_types = {"sel", "bkg", "gen", "sig", "eff", "ext", "dirt", "data", "mc_xsec", "data_xsec"};

    std::vector<std::string> vars = {"integrated", "reco_el_E", "true_el_E"
                                    //  "true_nu_E", "reco_nu_e"
                                     };
    
    std::vector<std::string> var_labels = {";;#nu_{e} + #bar{#nu}_{e} CC Cross-Section [10^{-39} cm^{2}]",
                                        ";Reco Leading Shower Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{reco}_{e}} CC Cross-Section [10^{-39} cm^{2}/GeV]",
                                        ";True Electron Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{true}_{e}} CC Cross-Section [10^{-39} cm^{2}/GeV]"
                                        // ";True #nu_{e} Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{true}_{#nu_{e}}} CC Cross-Section [10^{-39} cm^{2}/GeV]",
                                        // ";Reco #nu_{e} Energy [GeV];#frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{reco}_{#nu_{e}}} CC Cross-Section [10^{-39} cm^{2}/GeV]"
                                        };
    

    std::vector<std::string> reweighter_labels = {
        "CV",
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
        "weightsGenie",
        "weightsReint",
        "weightsPPFX"
    };


    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char *run_period, const char * xsec_file_in, Utility _utility, const char* run_mode);
    // -------------------------------------------------------------------------
    // Function to loop over events and calculate the cross section
    void LoopEvents(); 
    // -------------------------------------------------------------------------
    // Function to calculate the cross section
    double CalcCrossSec(double sel, double gen, double sig, double bkg, double flux, double ext, double dirt, double targ);
    // -------------------------------------------------------------------------
    // Function to calculate the cross section using binned histograms
    void CalcCrossSecHist(TH1D* h_sel, TH1D* h_eff, TH1D* h_bkg, double mc_scale_factor, double flux, double ext_scale_factor, TH1D* h_ext, double dirt_scale_factor, TH1D* h_dirt, TH1D* h_xsec, double targ, std::string mcdata);
    // -------------------------------------------------------------------------
    // Function to get the integrated flux
    double GetIntegratedFlux();
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
    // Initialise the histograms for this class
    void InitialiseHistograms(std::string run_mode);
    // -------------------------------------------------------------------------
    // Helper function to fill the histograms
    void FillHists(int label, int uni, int xsec_type, double weight_uni, float shr_energy_cali, float elec_e);
    // -------------------------------------------------------------------------
    // Function to return the integrated flux for universe i or the weight in case of a beamline variation
    double GetFluxUni(int uni, std::string value, std::string label);
    // -------------------------------------------------------------------------



    private:

    // Here we create the trees 



}; // End Class CrossSectionHelper

#endif