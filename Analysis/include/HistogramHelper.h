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
        k_reco_vtx_x,                                                           // Reco Vertex X
        k_reco_vtx_y,                                                           // Reco Vertex Y
        k_reco_vtx_z,                                                           // Reco Vertex Z
        k_reco_vtx_x_sce,                                                       // Reco Vertex X Space Charge Corrected
        k_reco_vtx_y_sce,                                                       // Reco Vertex Y Space Charge Corrected
        k_reco_vtx_z_sce,                                                       // Reco Vertex Z Space Charge Corrected
        k_reco_leading_mom,                                                     // Leading Shower Momentum
        k_reco_shower_to_vtx_dist,                                              // 2D distance shower vertex to reco nu vertex
        k_reco_track_to_vtx_dist,                                               // 2D distance track vertex to reco nu vertex
        k_reco_shr_hits_max,                                                    // Leading Shower hits in all planes
        k_reco_n_track_contained,                                               // Number of Tracks Contained
        k_reco_n_shower_contained,                                              // Number of Showers Contained
        k_reco_leading_shower_phi,                                              // Leading shower phi
        k_reco_shower_gamma,                                                    // Leading shower NuMI phi
        k_reco_leading_shower_theta,                                            // Leading shower theta
        k_reco_leading_shower_cos_theta,                                        // Leading shower cos theta
        k_reco_shower_multiplicity,                                             // Leading shower multiplicity
        k_reco_track_multiplicity,                                              // Leading track multiplicity
        k_reco_topological_score,                                               // Pandora Topological Score
        k_reco_track_shower_dist,                                               // Track shower dist
        k_reco_track_shower_angle,                                              // Track shower angle
        k_reco_hits_ratio,                                                      // Ratio hits from showers to slice
        k_reco_hits_ratio_th,                                                   // Ratio hits from showers to slice with a threshold
        k_reco_hits_ratio_ldg,                                                  // Ratio hits from leading shower to slice
        k_reco_hits_ratio_ldg_th,                                               // Ratio hits from leading shower to slice with a minimum hit requirement
        k_reco_shower_score,                                                    // Shower score
        k_reco_track_score,                                                     // Track score
        k_reco_effective_angle,                                                 // Angle between vector from target to vertex and shower direction.
        k_reco_effective_angle_rebin,                                           // Angle between vector from target to vertex and shower direction.
        k_reco_effective_cosangle,                                              // Cosine of the Angle between vector from target to vertex and shower direction.
        k_reco_effective_cosangle_rebin,                                        // Cosine of the Angle between vector from target to vertex and shower direction.
        k_reco_shower_energy_cali,                                              // Calibrated energy of just the leading shower 
        k_reco_shower_energy_cali_rebin,                                        // Calibrated energy of just the leading shower with optimised bins
        k_reco_shr_hits_tot,                                                    // Total number of hits for all showers
        k_reco_shr_hits_y_tot,                                                  // Total number of hits for all showers in the collection plane
        k_reco_shr_trkfit_2cm_dEdx_u,                                           // dE/dx of the leading shower on the Y plane with the track fitting, use first 2 cm u plane
        k_reco_shr_trkfit_2cm_dEdx_v,                                           // dE/dx of the leading shower on the Y plane with the track fitting, use first 2 cm v plane
        k_reco_shr_trkfit_2cm_dEdx_y,                                           // dE/dx of the leading shower on the Y plane with the track fitting, use first 2 cm y plane
        k_reco_shr_trkfit_gap05_dEdx_u,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 10mm u plane
        k_reco_shr_trkfit_gap05_dEdx_v,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 10mm v plane
        k_reco_shr_trkfit_gap05_dEdx_y,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 5mm y plane
        k_reco_shr_trkfit_gap10_dEdx_u,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 5mm u plane
        k_reco_shr_trkfit_gap10_dEdx_v,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 5mm v plane
        k_reco_shr_trkfit_gap10_dEdx_y,                                         // dE/dx of the leading shower on the Y plane with the track fitting, use 1x4 cm box, skip first 10mm y plane
        k_reco_opfilter_beam,                                                   // Common optical filter beam
        k_reco_opfilter_veto,                                                   // Common Optical filter michel veto
        k_reco_softwaretrig,                                                    // Software Trigger
        k_reco_nslice,                                                          // Pandora Slice ID
        k_reco_slclustfrac,                                                     // Reco Fraction of hits in the slice that are fully reconstructed to 3D particles.
        k_reco_cosmicIP,                                                        // Reco Closest distance between shower start and space points associated to tracks flagged as cosmics.
        k_reco_CosmicIPAll3D,                                                   // Reco 3D distance of shower start from closest spacepoint of any pfp not in the neutrino slice
        k_reco_CosmicDirAll3D,                                                  // cosine of 3D direction difference between shower and closest pfp not in the neutrino slice
        k_reco_shr_tkfit_dedx_u,                                                // The dEdx using the trackfit variable u plane
        k_reco_shr_tkfit_dedx_v,                                                // The dEdx using the trackfit variable v plane
        k_reco_shr_tkfit_dedx_y,                                                // The dEdx using the trackfit variable collection
        k_reco_shr_tkfit_dedx_max,                                              // The dEdx using the trackfit variable plane with the max hits
        k_reco_shr_tkfit_dedx_max_tune,                                         // The dEdx using the trackfit variable plane with the max hits with pi0 tune
        k_reco_shr_tkfit_dedx_max_with_tracks,                                  // The dEdx using the trackfit variable plane with the max hits for events with tracks only
        k_reco_shr_tkfit_dedx_y_no_tracks,                                      // The dEdx using the trackfit variable collection in the case there is no tracks
        k_reco_shr_tkfit_dedx_max_no_tracks,                                    // The dEdx using the trackfit variable plane with the most hits in the case there is no tracks
        k_reco_shr_tkfit_dedx_y_good_theta,                                     // The dEdx using the trackfit variable collection for angles not close to parallel to the y plane
        k_reco_shr_tkfit_dedx_y_bad_theta,                                      // The dEdx using the trackfit variable collection for angles close to parallel to the y plane
        k_reco_shr_tkfit_dedx_v_bad_theta,                                      // The dEdx using the trackfit variable v plane for angles close to parallel to the y plane
        k_reco_shr_tkfit_dedx_u_bad_theta,                                      // The dEdx using the trackfit variable v plane for angles close to parallel to the y plane
        k_reco_flash_time,                                                      // The Flash time
        k_reco_flash_pe,                                                        // The Flash PE
        k_reco_shrsubclusters,                                                  // Number of subclusters the shower can be broken into, Sum all three planes
        k_reco_shrmoliereavg,                                                   // Average angle between the shower’s direction and its 3D spacepoints.
        k_reco_shrmoliererms,                                                   // RMS of the moliere angle
        k_reco_shrPCA1CMed_5cm,                                                 // Median PCA component calculated in 5 cm blocks.
        k_reco_shrMCSMom,                                                       // Multiple Coulomb scattering shower momentum
        k_reco_closestNuCosmicDist,                                             // Distance between the neutrino vertex and (closest?) cosmic trajectory tagged from CRT
        k_reco_trk_len,                                                         // Length of the longest track
        k_reco_nu_e,                                                            // Reconstructed Neutrino Energy
        k_reco_contained_fraction,                                              // Ratio of PFP hits in FV to the slice
        k_reco_run_number,                                                      // Wont be used for stack, but for run normalisation plot
        k_reco_nu_purity_from_pfp,                                              // Purity
        k_reco_crtveto,                                                         // CRT veto
        k_reco_crthitpe,                                                        // CRT hit pe
        k_reco_shr_ang_numi,                                                    // Angle of the reconstructed leading shower relative to the numi beamline
        k_reco_single_bin,                                                      // Just fill with weight. We need this to get the uncertainty on the purity
        k_reco_trk_pid_score,                                                   // The PID score for tracks in the event
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
        k_eff_cosine_beta_nue,       // cosine beta nue
        k_eff_cosine_beta_nuebar,    // cosine beta nue bar
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
        k_true_elec_gamma,    // True phi of electron in NuMI coordinates
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
        k_true_elec_E_true_beta_rebin,
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