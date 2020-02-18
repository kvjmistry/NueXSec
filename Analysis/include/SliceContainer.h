#ifndef SLICECONTAINER_H
#define SLICECONTAINER_H

#include <vector>
#include <string>
#include <TTree.h>
#include <iostream>

#include "utility.h"

/* 
Class to hold information for the eventslice in the pandora ntuples for ease of use

*/

// Slice Container Class
class SliceContainer {
public:

    // -------------------------------------------------------------------------
    // Initialise the class
    void Initialise(TTree* tree);
    // -------------------------------------------------------------------------
    // Function to classify the slice
    std::pair<std::string, int>  SliceClassifier(int type);
    // -------------------------------------------------------------------------
    // Returns the Category defined by the pandora team 
    // (can be used as a cross check or making similar plots to them)
    std::string SliceCategory();
    // -------------------------------------------------------------------------
    // Function to return the genie interaction mode, e.g. ccqe, ccmec etc.
    std::string SliceInteractionType(int type);
    // -------------------------------------------------------------------------
    utility _util;

    int   selected;
    int   run;                   // Run
    int   sub;                   // Subrun
    int   evt;                   // Event
   
    // Shower Properties 
    float shr_energy_tot;        // Total Shower Energy
    float shr_energy;            // Shower Energy
    float shr_energy_tot_cali;   // Shower Sum of the energy of the calibrated showers (in GeV). Used only at pre-selection as a “Michel veto”.
    float shr_energy_cali;       // Shower Energy Cali
    float shr_theta;             // Shower Theta
    float shr_phi;               // Shower Phi
    float shr_pca_0;             // Shower PCA Plane 0
    float shr_pca_1;             // Shower PCA Plane 1
    float shr_pca_2;             // Shower PCA Plane 2
    float shr_px;                // Shower Px
    float shr_py;                // Shower Py
    float shr_pz;                // Shower Pz
    float shr_openangle;         // Shower Open Angle
    float shr_tkfit_start_x;     // Shower Track Fit Start X
    float shr_tkfit_start_y;     // Shower Track Fit Start Y
    float shr_tkfit_start_z;     // Shower Track Fit Start Z
    float shr_tkfit_theta;       // Shower Track Fit Theta
    float shr_tkfit_phi;         // Shower Track Fit Phi
    float shr_start_x;           // Shower Start X
    float shr_start_y;           // Shower Start Y
    float shr_start_z;           // Shower Start Z
    float shr_dedx_Y;            // Shower dEdx Y Plane 0
    float shr_dedx_V;            // Shower dEdx V Plane 1
    float shr_dedx_U;            // Shower dEdx U Plane 2
    float shr_dedx_Y_cali;       // Shower dEdx Cali Y Plane 0
    float shr_dedx_V_cali;       // Shower dEdx Cali V Plane 1
    float shr_dedx_U_cali;       // Shower dEdx Cali U Plane 2
    float shr_tkfit_dedx_Y;      // Shower TrackFit dEdx Cali Y Plane 0
    float shr_tkfit_dedx_V;      // Shower TrackFit dEdx Cali V Plane 1
    float shr_tkfit_dedx_U;      // Shower TrackFit dEdx Cali U Plane 2
    int   shr_tkfit_nhits_Y;     // Shower TrackFit Total Hits Y Plane 0
    int   shr_tkfit_nhits_V;     // Shower TrackFit Total Hits V Plane 1
    int   shr_tkfit_nhits_U;     // Shower TrackFit Total Hits U Plane 2
    float shr_tkfit_dedx_Y_alt;  // Shower TrackFit dEdx Y Plane 0 Alternate
    float shr_tkfit_dedx_V_alt;  // Shower TrackFit dEdx V Plane 1 Alternate
    float shr_tkfit_dedx_U_alt;  // Shower TrackFit dEdx U Plane 2 Alternate
    int   shr_tkfit_nhits_Y_alt; // Shower TrackFit Total Hits Y Plane 0 Alternate
    int   shr_tkfit_nhits_V_alt; // Shower TrackFit Total Hits V Plane 1 Alternate
    int   shr_tkfit_nhits_U_alt; // Shower TrackFit Total Hits U Plane 2 Alternate
    int   shr_tkfit_npoints;     // Shower TrackFit Number of Points
    float shr_trkfitmedangle;    // Shower TrackFit Median Angle
    float shrmoliereavg;         // Shower Average angle between the shower’s direction and its 3D spacepoints.
    float shrmoliererms;         // Shower Moliere RMS
    
    unsigned char  ismerged;     // Check if a proton is merged at the beginning of a shower.
    float merge_bestdot;
    float merge_bestdist; // Distance between shower start point and track start (or end) point for the track in the slice that best matches the direction of the shower.
    float merge_vtx_x;
    float merge_vtx_y;
    float merge_vtx_z;
    int   merge_tk_ipfp;
    
    float shr_tkfit_2cm_dedx_Y;   // Median dE/dx computed over the first 2 cm of the shower’s trunk. Plane Y
    float shr_tkfit_2cm_dedx_V;
    float shr_tkfit_2cm_dedx_U;
    unsigned int   shr_tkfit_2cm_nhits_Y;
    unsigned int   shr_tkfit_2cm_nhits_V;
    unsigned int   shr_tkfit_2cm_nhits_U;
    float shr_tkfit_gap05_dedx_Y;
    float shr_tkfit_gap05_dedx_V;
    float shr_tkfit_gap05_dedx_U;
    unsigned int   shr_tkfit_gap05_nhits_Y;
    unsigned int   shr_tkfit_gap05_nhits_V;
    unsigned int   shr_tkfit_gap05_nhits_U;
    float shr_tkfit_gap10_dedx_Y;  // Median dE/dx computed over the [1,5] cm of the shower’s trunk. Plane Y
    float shr_tkfit_gap10_dedx_V;
    float shr_tkfit_gap10_dedx_U;
    unsigned int   shr_tkfit_gap10_nhits_Y;
    unsigned int   shr_tkfit_gap10_nhits_V;
    unsigned int   shr_tkfit_gap10_nhits_U;
    
    float shr_chipr;            // Shower Chi Prob
    float shr_chimu;            // Shower Chi Mu
    float shr_bragg_p;          // Shower Bragg Likelihood P
    float shr_bragg_mu;         // Shower Bragg Likelihood Mu
    float shr_bragg_mip;        // Shower Bragg Likelihood MIP
    float shr_bragg_kaon;       // Shower Bragg Likelihood Kaon
    float shr_bragg_pion;       // Shower Bragg Likelihood Pion
    
    float tksh_distance;        // Distance between leading shower and longest track start points.
    float tksh_angle;           // Angle between leading shower and longest track directions.
    
    float shr_distance;         // Shower Distance between the shower start and the neutrino vertex. Labelled as shower_vtx_dist in technote
    float shr_score;            // Shower Pandora SVM track/shower score for the leading shower.
    int   shr_bkt_pdg;          // Shower Backtracked PDG
    float shr_bkt_purity;       // Shower Backtracked Purity
    float shr_bkt_completeness; // Shower Backtracked Completeness
    float shr_bkt_E;            // Shower Backtracked Energy
    
    float trk_len;              // Track Length
    float trk_theta;            // Track Theta
    float trk_phi;              // Track Phi 
    float trk_energy;           // Track Energy 
    float trk_energy_muon;      // Track Energy Muon
    float trk_energy_muon_mcs;  // Track Energy Muon MCS
    float trk_energy_tot;       // Track Energy Total
    float trk_energy_muon_tot;  // Track Energy Muon Total
    float trk_distance;         // Track Distance between the track start and the neutrino vertex.
    float trk_score;            // Track Score
    int   trk_bkt_pdg;          // Track Backtrack
    float trk_bkt_purity;       // Track Backtrack Purity
    float trk_bkt_completeness; // Track Backtrack Completeness
    float trk_bkt_E;            // Track Backtrack Energy
    float trk_chipr_best;       // Track Chi Prob Best
    float trk_chipr_worst;      // Track Chi Prob Worst
    float trk_chimu_best;       // Track Chi Mu Best
    float trk_chimu_worst;      // Track Chi Mu Worst
    float trk_chipr;            // Track Chi Probabilty
    float trk_chimu;            // Track Chi Mu
    float trk_pida;             // Track Chi PIDA
    
    float trk_bragg_p;         // Track Bragg Likelihood P
    float trk_bragg_mu;        // Track Bragg Likelihood Mu
    float trk_bragg_mip;       // Track Bragg Likelihood MIPP
    float trk_bragg_kaon;      // Track Bragg Likelihood Kaon
    float trk_bragg_pion;      // Track Bragg Likelihood Pion
    
    int   trk_hits_max;        // Track Maximum Hits
    int   shr_hits_max;        // Shower Maximum Hits
    
    float trkshrhitdist0;      // Track Shower Hit Distance Plane 0
    float trkshrhitdist1;      // Track Shower Hit Distance Plane 1
    float trkshrhitdist2;      // Track Shower Hit Distance Plane 2
    
    int   total_hits_y;        // Total Hits Y Plane 0
    float extra_energy_y;      // ????
    float trk_energy_hits_tot; // Track Energy Hits Total
    
    int   shrsubclusters0;
    int   shrsubclusters1;
    int   shrsubclusters2;
    float shrclusfrac0;
    float shrclusfrac1;
    float shrclusfrac2;
    float shrclusdir0;
    float shrclusdir1;
    float shrclusdir2;
    
    int   shr_hits_tot;         // Shower Hits Total
    int   shr_hits_y_tot;       // Shower Hits Total Y Plane 0
    int   shr_hits_u_tot;       // Shower Hits Total U Plane 1 
    int   shr_hits_v_tot;       // Shower Hits Total V Plane 2
    
    int   trk_hits_tot;         // Track Hits Total 
    int   trk_hits_y_tot;       // Track Hits Y Plane 0
    int   trk_hits_u_tot;       // Track Hits U Plane 2
    int   trk_hits_v_tot;       // Track Hits V Plane 1
    
    float elecclusters_U_charge;
    float elecclusters_V_charge;
    float elecclusters_Y_charge;
    int   elecclusters_U_N;
    int   elecclusters_V_N;
    int   elecclusters_Y_N;
    
    int   n_tracks_contained;    // Reco Number of tracks fully contained in the fiducial volume
    int   n_showers_contained;   // Reco Number of showers with a starting point within the fiducial volume.
    float matched_E;
    float hits_ratio;            // Reco Ratio between hits from showers and total number of hits in the slice.
    float contained_fraction;    // Reco Hits in PFParticles contained in the fiducial volume over the total number of clustered hits in the slice.
    float sps_contained_fraction;
    float pt;
    float p;
    float pt_assume_muon;
    float p_assume_muon;
    float dvtx;
    float dtrk;
    float contained_sps_ratio;
    
    float CosmicIP;          // Reco Closest distance between shower start and space points associated to tracks flagged as cosmics.
    float CosmicIPAll3D;
    float CosmicDirAll3D;
    float CosmicIPAll2DEnds;
    float CosmicDirAll2DEnds;
    float CosmicIPAll2DOvlp;
    float CosmicDirAll2DOvlp;
    
    float leeweight;
    
    float true_pt;
    float true_pt_visible;
    float true_p;
    float true_p_visible;
    float true_e_visible;
    
    float opfilter_pe_beam; // Common Optical Filter Decicion
    float opfilter_pe_veto;
    
    int   nu_pdg;      // True neutrino PDG
    int   ccnc;        // True CC or NC Interaction
    int   interaction; // True Interaction code from GENIE
    float nu_e;        // True Neutrino Energy [GeV]
    float nu_pt;       // True Neutrino Transverse Energy of Interaction
    float theta;       // True Neutrino Theta
    bool  isVtxInFiducial; // True Neutrino Vertex in FV
    
    // Is the truth information contained? 
    // Require all track start/end point in FV and showers deposit > 60% of energy
    // in TPC or deposit at least 100 MeV in TPC
    bool truthFiducial;
    
    float true_nu_vtx_t;         // True Neutrino Vtx t
    float true_nu_vtx_x;         // True Neutrino Vtx x
    float true_nu_vtx_y;         // True Neutrino Vtx y
    float true_nu_vtx_z;         // True Neutrino Vtx z
    float true_nu_vtx_sce_x;     // True Neutrino Vtx Space Charge x
    float true_nu_vtx_sce_y;     // True Neutrino Vtx Space Charge y
    float true_nu_vtx_sce_z;     // True Neutrino Vtx Space Charge z
    
    float reco_nu_vtx_x;         // Reco Neutrino Vtx x
    float reco_nu_vtx_y;         // Reco Neutrino Vtx y
    float reco_nu_vtx_z;         // Reco Neutrino Vtx z
    float reco_nu_vtx_sce_x;     // Reco Reconstructed neutrino interaction vertex in (x,y,z) coordinates. The space charged correction is applied
    float reco_nu_vtx_sce_y;     // Reco Neutrino Vtx Space Charge y
    float reco_nu_vtx_sce_z;     // Reco Neutrino Vtx Space Charge z
    
    int   nmuon;   // Number of Muons
    float muon_e;  // Muon Energy
    float muon_c;  // Muon Completeness
    float muon_p;  // Muon Purity
     
    int   nelec;   // Number of electrons
    float elec_e;  // Electron Energy
    float elec_c;  // Electron Completeness
    float elec_p;  // Electron Purity
    
    float elec_vx;   // Electron Vtx x
    float elec_vy;   // Electron Vtx y
    float elec_vz;   // Electron Vtx z
    float elec_px;   // Electron Px
    float elec_py;   // Electron Py
    float elec_pz;   // Electron Pz

    int   npi0;     // Number of Pi0
    float pi0_e;    // Pi0 Energy
    float pi0_c;    // Pi0 Completeness
    float pi0_p;    // Pi0 Purity
    
    int   nneutron; // Number of Neutrons
    
    int   nproton;  // Number of Protons
    float proton_e; // Proton Energy
    float proton_c; // Proton Completeness
    float proton_p; // Proton Purity
    
    int   npion;    // Number of Pions
    float pion_e;   // Pion Energy
    float pion_c;   // Pion Completeness
    float pion_p;   // Pion Purity
    
    int   nslice;     // Reco Number of neutrino slices identified by the SliceID. Values are 0 or 1.
    int   crtveto;    // Reco Boolean variable checking if the event passes the CRT veto.
    float crthitpe;   // CRT Hit PE
    int   category;   // Truth Category
    float lep_e;
    int   pass;
    int   swtrig;     // Software Trigger
    int   evnhits;
    int   slpdg;      // Reco Slice PDG of primary PFP
    int   slnhits;    // Reco Slice Number of Hits within it
    int   n_pfps;     // Reco Number of pfp's
    int   n_tracks;   // Reco Number of Tracks
    int   n_showers;  // Reco Number of Showers
    
    int   hits_u;
    int   hits_v;
    int   hits_y;
    
    float topological_score;  // Reco Topological Score
    float slclustfrac;        // Reco Fraction of hits in the slice that are fully reconstructed to 3D particles.
    float endmuonmichel;
    int   filter_antibdt;
    int   filter_ncpi0;
    int   filter_pi0;
    int   filter_ccinclusive;
    float flash_pe;                     // Reco Flash PE
    float flash_time;                   // Reco Flash Time
    float nu_flashmatch_score;          // Reco Neutrino Flashmatch Score
    float best_cosmic_flashmatch_score; // Reco Best Cosmic Flashmatch Score
    
    float NeutrinoEnergy0;
    float NeutrinoEnergy1;
    float NeutrinoEnergy2;
    
    float SliceCaloEnergy0;
    float SliceCaloEnergy1;
    float SliceCaloEnergy2;
    
    int   nflag_pl1;
    int   nnoise_pl1;
    int   nslhits_pl1;
    int   nslnoise_pl1;
    int   nhits_pl1;
    float frac_slnoise_pl1;
    
    float secondshower_U_charge;
    int   secondshower_U_nhit;
    float secondshower_U_vtxdist;
    float secondshower_U_eigenratio;
    float secondshower_U_dot;
    float secondshower_U_dir;
    float secondshower_V_charge;
    int   secondshower_V_nhit;
    float secondshower_V_vtxdist;
    float secondshower_V_eigenratio;
    float secondshower_V_dot;
    float secondshower_V_dir;
    float secondshower_Y_charge;
    int   secondshower_Y_nhit;
    float secondshower_Y_vtxdist;
    float secondshower_Y_eigenratio;
    float secondshower_Y_dot;
    float secondshower_Y_dir;
    
    int   evnunhits;
    int   evlepnhits;
    int   evpronhits;
    int   evpi1nhits;
    int   evpi0nhits;
    int   evneunhits;
    int   evgamnhits;
    int   evothnhits;
    
    int   slnunhits;   // Slice Neutrino Number of Hits
    int   sllepnhits;  // Slice Lepton Number of Hits
    int   slpronhits;  // Slice Proton Nummber of Hits
    int   slpi1nhits;  // Slice 1 Pion Number of Hits
    int   slpi0nhits;  // Slice 0 Pion Number of Hits
    int   slneunhits;  // Slice Neutrino Number of Hits
    int   slgamnhits;  // Slice Gamma Number of Hits
    int   slothnhits;  // Slice Other Number of Hits
    
    float nu_completeness_from_pfp; // Neutrino Completeness from PFP ()
    float nu_purity_from_pfp;       // Neutrino Purity from PFP (how many out of all the hits are the neutrino) (needs to be higher than 50% otherwise tag as a mixed)
    int   n_tracks_pandora;         // Number of Tracks Returned by Pandora
    
    float vtx_fit_pandora_x;
    float vtx_fit_pandora_y;
    float vtx_fit_pandora_z;
    
    int   n_tracks_tkfit;
    float vtx_fit_tkfit_x;
    float vtx_fit_tkfit_y;
    float vtx_fit_tkfit_z;
    
    float bdt_nuNCpi0;
    float bdt_numuCCpi0;
    float bdt_numuCC;
    float bdt_ext;
    float bdt_cosmic;
    float bdt_global;

    std::vector<float> *pfp_generation_v          = nullptr;
    std::vector<float> *pfp_trk_daughters_v       = nullptr;
    std::vector<float> *pfp_shr_daughters_v       = nullptr;
    std::vector<float> *trk_score_v               = nullptr;
    std::vector<float> *pfpdg_v                   = nullptr;
    std::vector<float> *pfnhits_v                 = nullptr;
    std::vector<float> *pfnplanehits_U_v          = nullptr;
    std::vector<float> *pfnplanehits_V_v          = nullptr;
    std::vector<float> *pfnplanehits_Y_v          = nullptr;
    std::vector<float> *pfpplanesubclusters_U_v   = nullptr;
    std::vector<float> *pfpplanesubclusters_V_v   = nullptr;
    std::vector<float> *pfpplanesubclusters_Y_v   = nullptr;
    std::vector<float> *pfpplanesubhitfracmax_U_v = nullptr;
    std::vector<float> *pfpplanesubhitfracmax_V_v = nullptr;
    std::vector<float> *pfpplanesubhitfracmax_Y_v = nullptr;
    
    std::vector<int>   *mc_pdg_v  = nullptr;
    std::vector<float> *mc_E_v    = nullptr;
    std::vector<float> *mc_vx_v   = nullptr;
    std::vector<float> *mc_vy_v   = nullptr;
    std::vector<float> *mc_vz_v   = nullptr;
    std::vector<float> *mc_endx_v = nullptr;
    std::vector<float> *mc_endy_v = nullptr;
    std::vector<float> *mc_endz_v = nullptr;
    std::vector<float> *mc_px_v   = nullptr;
    std::vector<float> *mc_py_v   = nullptr;
    std::vector<float> *mc_pz_v   = nullptr;
    std::vector<float> *mc_completeness_v = nullptr;
    std::vector<float> *mc_purity_v       = nullptr;
    
    std::map<std::string, std::vector<double>> *weights_v = nullptr;
    std::vector<double> *weightsFlux_v   = nullptr;
    std::vector<double> *weightsGenie_v  = nullptr;
    std::vector<double> *weightsReint_v  = nullptr;
    
    std::vector<float> *cosmic_flashmatch_score_v = nullptr;
    std::vector<float> *peSpectrum_v              = nullptr;
    std::vector<float> *peHypothesisNu_v          = nullptr;
    std::vector<float> *peHypothesisCosmic_v      = nullptr;
    
    std::vector<float> *shr_dedx_u_v        = nullptr;
    std::vector<float> *shr_dedx_v_v        = nullptr;
    std::vector<float> *shr_dedx_y_v        = nullptr;
    std::vector<float> *shr_energy_u_v      = nullptr;
    std::vector<float> *shr_energy_v_v      = nullptr;
    std::vector<float> *shr_energy_y_v      = nullptr;
    std::vector<size_t> *shr_pfp_id_v       = nullptr;
    std::vector<float> *shr_start_x_v       = nullptr;
    std::vector<float> *shr_start_y_v       = nullptr;
    std::vector<float> *shr_start_z_v       = nullptr;
    std::vector<float> *shr_dist_v          = nullptr;
    std::vector<float> *shr_start_U_v       = nullptr;
    std::vector<float> *shr_start_V_v       = nullptr;
    std::vector<float> *shr_px_v            = nullptr;
    std::vector<float> *shr_py_v            = nullptr;
    std::vector<float> *shr_pz_v            = nullptr;
    std::vector<float> *shr_openangle_v     = nullptr;
    std::vector<float> *shr_theta_v         = nullptr;
    std::vector<float> *shr_phi_v           = nullptr;
    std::vector<float> *shr_pitch_u_v       = nullptr;
    std::vector<float> *shr_pitch_v_v       = nullptr;
    std::vector<float> *shr_pitch_y_v       = nullptr;
    std::vector<int> *shr_tkfit_nhits_v     = nullptr;
    std::vector<float> *shr_tkfit_start_x_v = nullptr;
    std::vector<float> *shr_tkfit_start_y_v = nullptr;
    std::vector<float> *shr_tkfit_start_z_v = nullptr;
    std::vector<float> *shr_tkfit_start_U_v = nullptr;
    std::vector<float> *shr_tkfit_start_V_v = nullptr;
    std::vector<float> *shr_tkfit_theta_v   = nullptr;
    std::vector<float> *shr_tkfit_phi_v     = nullptr;
    std::vector<float> *shr_tkfit_pitch_u_v = nullptr;
    std::vector<float> *shr_tkfit_pitch_v_v = nullptr;
    std::vector<float> *shr_tkfit_pitch_y_v = nullptr;
    std::vector<float> *shr_tkfit_dedx_u_v  = nullptr;
    std::vector<float> *shr_tkfit_dedx_v_v  = nullptr;
    std::vector<float> *shr_tkfit_dedx_y_v  = nullptr;
    std::vector<float> *shr_tkfit_gap10_dedx_u_v = nullptr;
    std::vector<float> *shr_tkfit_gap10_dedx_v_v = nullptr;
    std::vector<float> *shr_tkfit_gap10_dedx_y_v = nullptr;
    std::vector<float> *shr_tkfit_dedx_nhits_u_v = nullptr;
    std::vector<float> *shr_tkfit_dedx_nhits_v_v = nullptr;
    std::vector<float> *shr_tkfit_dedx_nhits_y_v = nullptr;
    std::vector<float> *shr_llr_pid_u_v     = nullptr;
    std::vector<float> *shr_llr_pid_v_v     = nullptr;
    std::vector<float> *shr_llr_pid_y_v     = nullptr;
    std::vector<float> *shr_llr_pid_v       = nullptr;
    std::vector<float> *shr_llr_pid_score_v = nullptr;
    std::vector<float> *shr_moliere_avg_v   = nullptr;
    std::vector<float> *shr_moliere_rms_v   = nullptr;
    // std::vector<float> *shr_spacepoint_start_x_v = nullptr;
    // std::vector<float> *shr_spacepoint_start_y_v = nullptr;
    // std::vector<float> *shr_spacepoint_start_z_v = nullptr;
    // std::vector<float> *shr_spacepoint_start_U_v = nullptr;
    // std::vector<float> *shr_spacepoint_start_V_v = nullptr;
    // std::vector<float> *shr_hits_start_U_wire_v = nullptr;
    // std::vector<float> *shr_hits_start_U_x_v = nullptr;
    // std::vector<float> *shr_hits_start_V_wire_v = nullptr;
    // std::vector<float> *shr_hits_start_V_x_v = nullptr;
    // std::vector<float> *shr_hits_start_Y_wire_v = nullptr;
    // std::vector<float> *shr_hits_start_Y_x_v = nullptr;

    std::vector<float> *trk_bragg_p_v     = nullptr;
    std::vector<float> *trk_bragg_mu_v    = nullptr;
    std::vector<float> *trk_bragg_mip_v   = nullptr;
    std::vector<float> *trk_pida_v        = nullptr;
    std::vector<float> *trk_pid_chipr_v   = nullptr;
    std::vector<float> *trk_pid_chipi_v   = nullptr;
    std::vector<float> *trk_pid_chika_v   = nullptr;
    std::vector<float> *trk_pid_chimu_v   = nullptr;
    std::vector<float> *trk_bragg_p_u_v   = nullptr;
    std::vector<float> *trk_bragg_mu_u_v  = nullptr;
    std::vector<float> *trk_bragg_mip_u_v = nullptr;
    std::vector<float> *trk_pida_u_v      = nullptr;
    std::vector<float> *trk_pid_chipr_u_v = nullptr;
    std::vector<float> *trk_pid_chipi_u_v = nullptr;
    std::vector<float> *trk_pid_chika_u_v = nullptr;
    std::vector<float> *trk_pid_chimu_u_v = nullptr;
    std::vector<float> *trk_bragg_p_v_v   = nullptr;
    std::vector<float> *trk_bragg_mu_v_v  = nullptr;
    std::vector<float> *trk_bragg_mip_v_v = nullptr;
    std::vector<float> *trk_pida_v_v      = nullptr;
    std::vector<float> *trk_pid_chipr_v_v = nullptr;
    std::vector<float> *trk_pid_chipi_v_v = nullptr;
    std::vector<float> *trk_pid_chika_v_v = nullptr;
    std::vector<float> *trk_pid_chimu_v_v = nullptr;
    std::vector<size_t> *trk_pfp_id_v     = nullptr;
    std::vector<float> *trk_dir_x_v       = nullptr;
    std::vector<float> *trk_dir_y_v       = nullptr;
    std::vector<float> *trk_dir_z_v       = nullptr;
    std::vector<float> *trk_start_x_v     = nullptr;
    std::vector<float> *trk_start_y_v     = nullptr;
    std::vector<float> *trk_start_z_v     = nullptr;
    std::vector<float> *trk_sce_start_x_v = nullptr;
    std::vector<float> *trk_sce_start_y_v = nullptr;
    std::vector<float> *trk_sce_start_z_v = nullptr;
    std::vector<float> *trk_end_x_v       = nullptr;
    std::vector<float> *trk_end_y_v       = nullptr;
    std::vector<float> *trk_end_z_v       = nullptr;
    std::vector<float> *trk_sce_end_x_v   = nullptr;
    std::vector<float> *trk_sce_end_y_v   = nullptr;
    std::vector<float> *trk_sce_end_z_v   = nullptr;
    std::vector<float> *trk_distance_v    = nullptr;
    std::vector<float> *trk_theta_v       = nullptr;
    std::vector<float> *trk_phi_v         = nullptr;
    std::vector<float> *trk_len_v         = nullptr;
    std::vector<float> *trk_mcs_muon_mom_v   = nullptr;
    std::vector<float> *trk_range_muon_mom_v = nullptr;
    std::vector<float> *trk_energy_proton_v  = nullptr;
    std::vector<float> *trk_energy_muon_v    = nullptr;
    std::vector<float> *trk_calo_energy_u_v  = nullptr;
    std::vector<float> *trk_calo_energy_v_v  = nullptr;
    std::vector<float> *trk_calo_energy_y_v  = nullptr;
    std::vector<float> *trk_llr_pid_u_v      = nullptr;
    std::vector<float> *trk_llr_pid_v_v      = nullptr;
    std::vector<float> *trk_llr_pid_y_v      = nullptr;
    std::vector<float> *trk_llr_pid_v        = nullptr;
    std::vector<float> *trk_llr_pid_score_v  = nullptr;
    
};

#endif