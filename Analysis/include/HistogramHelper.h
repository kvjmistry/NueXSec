#ifndef HISTOGRAMHELPER_H
#define HISTOGRAMHELPER_H

#include "SelectionCuts.h"
#include "SliceContainer.h"

// Class for filling and saving histograms. 
class HistogramHelper{

    public:
    // Default constructor
    HistogramHelper(){};
    
    // Destructor 
    ~HistogramHelper(); 

    // The output file
    TFile* f_nuexsec;

    // Class instances
    Utility _util;
    SelectionCuts _scuts;
    int _type{1};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type, const char * file_out, Utility util );
    // -------------------------------------------------------------------------
    // Initialise histograms
    void InitHistograms();
    // -------------------------------------------------------------------------
    // Function to make the directory hirarchy
    void MakeDirectory();
    // -------------------------------------------------------------------------
    // Function to fill the reco variables
    void FillHists(int type, int classification_index, std::string interaction, int _par_type, int cut_index, SliceContainer SC, double weight);
    // -------------------------------------------------------------------------
    // Function to write the histograms to a file
    void WriteReco(int type);
    // -------------------------------------------------------------------------
    // As above, but for histograms by particle type
    void WriteRecoPar(int type);
    // -------------------------------------------------------------------------
    // Function to fill the true neutrino in FV graphs
    void FillTEfficiency(int cut_index, std::string classification, SliceContainer SC, double weight);
    // -------------------------------------------------------------------------
    // Function to write the TEfficiency Graphs to file
    void WriteTEfficiency();
    // -------------------------------------------------------------------------
    // Function to write the truth histograms to file
    void WriteTrue();
    // -------------------------------------------------------------------------
    // Function to write the flash histograms
    void WriteFlash();
    // -------------------------------------------------------------------------
    // Function to write the interaction histograms
    void WriteInteractions();
    // -------------------------------------------------------------------------
    // Function to write the 2D signal vs Background Histograms
    void Write_2DSigBkgHists();
    // -------------------------------------------------------------------------
    // Pi Zero Stacked Histogram
    void FillPiZeroHists(int classification_index, SliceContainer SC, double weight, int pizero_mode);
    // -------------------------------------------------------------------------
    // Write the PiZero Histograms
    void WritePiZero(int type);
    // -------------------------------------------------------------------------
    // Write the NuMi Histograms
    void WriteNuMu(int type);
    // -------------------------------------------------------------------------
    // NuMu Stacked Histogram
    void FillNuMuHists(int classification_index, SliceContainer SC, double weight);
    // -------------------------------------------------------------------------



    private:

    // Here we create the histograms
   
    // vector of histograms to make, indexed by enums
    std::vector<std::vector<std::vector<TH1D*>>> TH1D_hists; 

    // vector of 2D histograms to make, indexed by enums
    std::vector<std::vector<std::vector<TH2D*>>> TH2D_hists_cuts; 

    // vector of histograms to make, indexed by enums -- for particle type
    std::vector<std::vector<std::vector<TH1D*>>> TH1D_hists_particle; 

    // Histograms for pi0 
    std::vector<std::vector<TH1D*>> TH1D_pi0_hists;

    // Histograms for NuMu
    std::vector<std::vector<TH1D*>> TH1D_numu_hists;

    // Histograms for the efficiency plot
    std::vector<std::vector<TH1D*>> TEfficiency_hists;

    // True histograms
    std::vector<std::vector<TH1D*>> TH1D_true_hists;
    std::vector<std::vector<TH2D*>> TH2D_true_hists;
    
    // Flash Histograms
    std::vector<TH1D*> TH1D_flash_hists;

    // Interaction Histograms
    std::vector<std::vector<std::vector<TH1D*>>> TH1D_interaction_hists; // sum/nue/nubar -- unselected/selected -- interaction type

    // 2D histograms for Signal and Background Rejection
    std::vector<std::vector<TH2D*>> TH2D_hists;
    
    // enum for histogram vars
    enum TH1D_hist_vars {
        k_reco_pi0mass,                                                         // What does the pi0 mass peak look like in a nuecc enriched sample
        k_TH1D_MAX
    };

    // enum for histogram vars
    enum TH1D_par_hist_vars {
        k_reco_shr_tkfit_dedx_max_par,  // dEdx in the plane with most hits with trackfit 1x4 cm box
        k_reco_shr_tkfit_dedx_y_par,    // dEdx in the collection plane with trackfit 1x4 cm box
        k_reco_trk_pid_score_par,           // Track PID Score
        k_TH1D_par_MAX
    };

    // enum for plots by efficiency
    enum TH1D_eff_vars {
        k_eff_nu_E,                  // True Electron-neutrino energy
        k_eff_elec_E,                // True Electron Energy
        k_eff_elec_E_many_bins,      // True Electron Energy with many bins
        k_eff_elec_E_rebin,          // True energy of electron with binning scheme
        k_eff_elec_E_rebin_nue,      // True energy of electron with binning scheme (nue only)
        k_eff_elec_E_rebin_nuebar,   // True energy of positron with binning scheme (nuebar only)
        k_eff_nu_E_nue,              // True Electron-neutrino energy
        k_eff_nu_E_nuebar,           // True anti Electron-neutrino energy
        k_eff_nu_E_single_bin,       // True Electron-neutrino energy, single bin
        k_eff_nu_E_nue_single_bin,   // True Electron-neutrino energy single bin
        k_eff_nu_E_nuebar_single_bin,// True anti Electron-neutrino energy single bin
        k_eff_nu_flash_time,         // Efficiency as a function of flash time
        k_eff_nu_theta,              // Efficiency as a function of nu theta
        k_eff_nu_phi,                // Efficiency as a function of nu phi
        k_eff_elec_theta,            // Efficiency as a function of electron theta
        k_eff_elec_phi,              // Efficiency as a function of electron phi
        k_eff_beta,                  // beta
        k_eff_beta_rebin,            // beta with uneven bins
        k_eff_beta_rebin_nue,        // beta with uneven bins nue
        k_eff_beta_rebin_nuebar,     // beta with uneven bins nueber
        k_eff_cosine_beta,           // cosine beta
        k_eff_cosine_beta_rebin,     // cosine beta with uneven bins
        k_eff_cosine_beta_rebin_nue, // cosine beta with uneven bins
        k_eff_cosine_beta_rebin_nuebar, // cosine beta with uneven bins
        k_eff_proton_multi,          // Efficiency as a function of the number of true protons in the interaction.
        k_eff_proton_multi_nue,      // Efficiency as a function of the number of true protons in the interaction. Nue events only
        k_eff_proton_multi_nuebar,   // Efficiency as a function of the number of true protons in the interaction. Nuebar events only
        k_eff_pion_multi,            // Efficiency as a function of the number of true pions in the interaction.
        k_eff_pion_multi_nue,        // Efficiency as a function of the number of true pions in the interaction. Nue events only
        k_eff_pion_multi_nuebar,     // Efficiency as a function of the number of true pions in the interaction. Nuebar events only
        k_eff_charg_par_multi,       // Efficiency as a function of the number of true charged particles (protons + pions) in the interaction.
        k_eff_charg_par_multi_nue,   // Efficiency as a function of the number of true charged particles (protons + pions) in the interaction. Nue events only
        k_eff_charg_par_multi_nuebar,// Efficiency as a function of the number of true charged particles (protons + pions) in the interaction. Nuebar events only
        k_TH1D_eff_MAX
    };

    enum TH1D_true_hist_vars {
        k_true_nue_theta,     // True nue theta in BNB coordinates (up from beam dir)
        k_true_nue_phi,       // True nue phi in BNB coordinates (around beam dir)
        k_true_nue_angle,     // True nue angle from numi beamline 
        k_true_nue_px,        // True nue px
        k_true_nue_py,        // True nue py
        k_true_nue_pz,        // True nue pz
        k_true_nue_e,         // True nue energy
        k_true_nue_p,         // True nue momentum
        k_true_vtx_x,         // True Vertex X
        k_true_vtx_y,         // True Vertex Y
        k_true_vtx_z,         // True Vertex Z
        k_true_vtx_x_sce,     // True Vertex X Space Charge Corrected
        k_true_vtx_y_sce,     // True Vertex Y Space Charge Corrected
        k_true_vtx_z_sce,     // True Vertex Z Space Charge Corrected
        k_true_elec_ang_targ, // True angle of electron shower wrt target
        k_true_elec_E,        // True energy of electron
        k_true_elec_theta,    // True theta of electron in BNB coordinates
        k_true_elec_phi,      // True phi of electron in BNB coordinates
        k_true_nu_ang_targ,   // True angle of electron shower wrt target
        k_reco_true_ang,      // Angle between the reco and true neutrino angle
        k_TH1D_true_MAX
    };

    enum TH1D_flash_hist_vars {
        k_flash_time,
        k_flash_time_single_bin,
        k_flash_pe,
        k_flash_time_sid1,   // Slice Id 1 (neutrino candidates)
        k_flash_pe_sid1,
        k_flash_time_sid0,  // Slice Id 0 (non neutrino candidates)
        k_flash_pe_sid0,
        k_TH1D_flash_MAX
    };

    enum TH2D_true_hist_vars {
        k_true_nue_phi_theta, 
        k_true_nue_energy_theta,
        k_true_nue_energy_phi,
        k_true_nue_energy_angle,
        k_true_nue_vtx_z_y,
        k_true_nue_vtx_z_y_sce,
        k_true_elec_E_reco_elec_E,
        k_true_nu_E_reco_nu_E,
        k_true_elec_E_reco_elec_E_extra_bins,
        k_true_elec_E_reco_elec_E_extra_bins_nue,
        k_true_elec_E_reco_elec_E_extra_bins_nuebar,
        k_true_nu_E_reco_nu_E_extra_bins,
        k_true_nu_vtx_x_reco_nu_vtx_x,
        k_true_nu_vtx_y_reco_nu_vtx_y,
        k_true_nu_vtx_z_reco_nu_vtx_z,
        k_true_shr_energy_purity,           // Actually purity as a function of reco shower
        k_true_shr_energy_completeness,     // Actually completeness as a function of reco shower
        k_true_shr_energy_resolution_reco,  // Actually resolution normed to reco as a function of reco shower
        k_true_shr_energy_resolution_true,  // Actually resolution normed to true as a function of reco shower
        k_true_shr_cosbeta_purity,           // Actually purity as a function of reco shower
        k_true_shr_cosbeta_completeness,     // Actually completeness as a function of reco shower
        k_true_shr_cosbeta_resolution_reco,  // Actually resolution normed to reco as a function of reco shower
        k_true_shr_cosbeta_resolution_true,  // Actually resolution normed to true as a function of reco shower
        k_elec_true_beta_reco_beta,
        k_elec_true_beta_reco_beta_nue,
        k_elec_true_beta_reco_beta_nuebar,
        k_elec_true_theta_reco_theta,
        k_elec_true_phi_reco_phi,
        k_elec_true_cosbeta_reco_cosbeta_rebin, 
        k_true_elec_E_reco_elec_E_rebin,
        k_TH2D_true_MAX
    };

    // 2D Histograms for separating signal and background
    enum TH2D_reco_hist_vars {
        k_reco_shr_dEdx_shr_dist,            // dedx y vs shr vtx distance
        k_reco_shr_dEdx_shr_dist_post,       // after the cut
        k_reco_shr_dEdx_max_shr_dist,        // Using max variable rather than just collection plane
        k_reco_shr_dEdx_max_shr_dist_post,   // after the cut
        k_reco_shr_dEdx_shr_dist_large_dedx, // for dedx values > 10 MeV/cm
        k_reco_shr_dEdx_moliere,             // dedx y and moliere average
        k_reco_shr_moliere_shr_dist,         // moliere average and shr vertex distance
        k_TH2D_reco_MAX
    };

    enum TH1D_pi0_hist_vars {
        k_pi0_mass,             // The pi0 mass peak no weighting 
        k_pi0_mass_norm,        // The pi0 mass peak normalisation fix
        k_pi0_mass_EScale,      // The pi0 mass peak energy dependent scaling
        k_pi0_energy,           // The pi0 energy dist no weighting 
        k_pi0_energy_norm,      // The pi0 energy dist normalisation fix
        k_pi0_energy_EScale,    // The pi0 energy dist energy dependent scaling
        k_TH1D_pi0_MAX
    };

    enum TH1D_numu_hist_vars {
        k_track_theta,      // Longest track theta 
        k_track_cos_theta,  // Longest track cos theta
        k_track_phi,        // Longest track phi
        k_muon_topo_score,  // Topological score (after muon selection)
        k_TH1D_numu_MAX
    };

    // Define the enums for 2D histograms broken down by cuts and classifications
    enum TH2D_cut_vars {
        k_2D_dedx_shower_energy,      // 2D plot of dedx and reconstrcted shower energy
        k_TH2D_cut_MAX
    };

    enum TH1D_interaction_vars {
        k_int_nu_E_nue,
        k_int_nu_E_nuebar,
        k_int_nu_E_nue_nuebar,
        k_int_nu_E_single_bin,
        k_int_nu_E_nue_single_bin,
        k_int_nu_E_nuebar_single_bin,
        k_int_elec_E,
        k_int_elec_E_nue,
        k_int_elec_E_nuebar,
        k_int_elec_E_rebin,
        k_int_elec_E_rebin_nue,
        k_int_elec_E_rebin_nuebar,
        k_int_elec_theta,
        k_int_elec_phi,
        k_int_effective_ang,
        k_int_beta_nue,
        k_int_beta_nuebar,
        k_int_cosbeta,
        k_int_cosbeta_rebin_nue,
        k_int_cosbeta_rebin_nuebar,
        k_INTERACTION_MAX
    };

}; // End Class Histogram Helper 

#endif