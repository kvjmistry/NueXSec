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
    void Initialise(TTree* tree, int type, TFile *f_flux_weights, const char * _run_period, utility util);
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
    // Function to classify the event by particle type of the leading shower
    std::pair<std::string, int> ParticleClassifier(int type);
    // -------------------------------------------------------------------------
    // Function to Get the PPFX CV correction weight
    double GetPPFXCVWeight();
    // -------------------------------------------------------------------------

    utility _util;

    int   selected;              // Adds at least 1 shower and containment????
    int   run;                   // Run
    int   sub;                   // Subrun
    int   evt;                   // Event
    std::string run_period;      // Run period 
    
    // Shower Properties 
    float shr_energy_tot;        // Shower: the energy of the showers (in GeV)
    float shr_energy;            // Shower: Energy of the shower with the largest number of hits (in GeV)
    float shr_energy_tot_cali;   // Shower: Sum of the energy of the calibrated showers (in GeV). Used only at pre-selection as a “Michel veto”.  does not take into account pi0. Only accounts for the fact the shower could be broken up.
    float shr_energy_cali;       // Shower: Energy of the calibrated shower with the largest number of hits (in GeV)
    
    float shr_theta;             // Shower: Reconstructed theta angle for the leading shower
    float shr_phi;               // Shower: Reconstructed phi angle for the leading shower
    
    float shr_pca_0;             // Shower: First eigenvalue of the PCAxis of the leading shower
    float shr_pca_1;             // Shower: Second eigenvalue of the PCAxis of the leading shower
    float shr_pca_2;             // Shower: Third eigenvalue of the PCAxis of the leading shower
    
    float shr_px;                // Shower: X component of the reconstructed momentum of the leading shower (in GeV/c)
    float shr_py;                // Shower: Y component of the reconstructed momentum of the leading shower (in GeV/c)
    float shr_pz;                // Shower: Z component of the reconstructed momentum of the leading shower (in GeV/c)
    
    float shr_openangle;         // Shower: Opening angle of the shower -- variable does not work...
    
    float shr_tkfit_start_x;     // Shower: Start x coordinate of the leading shower obtained with the track fitting
    float shr_tkfit_start_y;     // Shower: Start y coordinate of the leading shower obtained with the track fitting
    float shr_tkfit_start_z;     // Shower: Start z coordinate of the leading shower obtained with the track fitting
    
    float shr_tkfit_theta;       // Shower: Track angle of the leading shower obtained with the track fitting
    float shr_tkfit_phi;         // Shower: Phi angle of the leading shower obtained with the track fitting
    
    float shr_start_x;           // Shower: Start x coordinate of the leading shower
    float shr_start_y;           // Shower: Start y coordinate of the leading shower
    float shr_start_z;           // Shower: Start z coordinate of the leading shower
    
    float shr_dedx_Y;            // Shower: dE/dx of the leading shower on the Y plane with the 1x4 cm box method| Plane 2
    float shr_dedx_V;            // Shower: dE/dx of the leading shower on the V plane with the 1x4 cm box method| Plane 1
    float shr_dedx_U;            // Shower: dE/dx of the leading shower on the U plane with the 1x4 cm box method| Plane 0
    
    float shr_dedx_Y_cali;       // Shower: Calibrated dE/dx of the leading shower on the Y plane with the 1x4 cm box method| Plane 2
    float shr_dedx_V_cali;       // Shower: Calibrated dE/dx of the leading shower on the V plane with the 1x4 cm box method| Plane 1
    float shr_dedx_U_cali;       // Shower: Calibrated dE/dx of the leading shower on the U plane with the 1x4 cm box method| Plane 0
    
    float shr_tkfit_dedx_Y;      // Shower: dE/dx of the leading shower on the Y plane with the track fitting
    float shr_tkfit_dedx_V;      // Shower: dE/dx of the leading shower on the V plane with the track fitting
    float shr_tkfit_dedx_U;      // Shower: dE/dx of the leading shower on the U plane with the track fitting
    
    int   shr_tkfit_nhits_Y;     // Shower: Number of hits in the 1x4 cm box on the Y plane with the track fitting
    int   shr_tkfit_nhits_V;     // Shower: Number of hits in the 1x4 cm box on the V plane with the track fitting
    int   shr_tkfit_nhits_U;     // Shower: Number of hits in the 1x4 cm box on the U plane with the track fitting
    
    float shr_tkfit_dedx_Y_alt;  // Shower: dE/dx of the leading shower on the Y plane with the track fitting  [calculated using XYZ instead of Residual Range]
    float shr_tkfit_dedx_V_alt;  // Shower: dE/dx of the leading shower on the V plane with the track fitting  [calculated using XYZ instead of RR]
    float shr_tkfit_dedx_U_alt;  // Shower: dE/dx of the leading shower on the U plane with the track fitting  [calculated using XYZ instead of RR]
    
    int   shr_tkfit_nhits_Y_alt; // Shower: Number of hits in the 1x4 cm box on the Y plane with the track fitting [calculated using XYZ instead of RR]
    int   shr_tkfit_nhits_V_alt; // Shower: Number of hits in the 1x4 cm box on the V plane with the track fitting [calculated using XYZ instead of RR]
    int   shr_tkfit_nhits_U_alt; // Shower: Number of hits in the 1x4 cm box on the U plane with the track fitting [calculated using XYZ instead of RR]
    
    int   shr_tkfit_npoints;     // Shower: TrackFit Number of Points
    float shr_trkfitmedangle;    // Shower: TrackFit Median Angle
    float shrmoliereavg;         // Shower: Average angle between the shower’s direction and its 3D spacepoints.
    float shrmoliererms;         // Shower: rms of moliere angle
    
    unsigned char  ismerged;     // Check if a proton is merged at the beginning of a shower.
    float merge_bestdot;
    float merge_bestdist; // Distance between shower start point and track start (or end) point for the track in the slice that best matches the direction of the shower.
    float merge_vtx_x;
    float merge_vtx_y;
    float merge_vtx_z;
    int   merge_tk_ipfp;
    
    float shr_tkfit_2cm_dedx_Y;   // dE/dx of the leading shower on the Y plane with the track fitting, use first 2 cm
    float shr_tkfit_2cm_dedx_V;
    float shr_tkfit_2cm_dedx_U;
    
    unsigned int   shr_tkfit_2cm_nhits_Y; // Number of hits in the 1x4 cm box on the Y plane with the track fitting, use first 2 cm
    unsigned int   shr_tkfit_2cm_nhits_V;
    unsigned int   shr_tkfit_2cm_nhits_U;
    
    float shr_tkfit_gap05_dedx_Y; // dE/dx of the leading shower on the Y plane with the track fitting, skip first 5 mm
    float shr_tkfit_gap05_dedx_V;
    float shr_tkfit_gap05_dedx_U;
    
    unsigned int   shr_tkfit_gap05_nhits_Y; //Number of hits in the 1x4 cm box on the Y plane with the track fitting, skip first 5 mm
    unsigned int   shr_tkfit_gap05_nhits_V;
    unsigned int   shr_tkfit_gap05_nhits_U;
    
    float shr_tkfit_gap10_dedx_Y;  // dE/dx of the leading shower on the Y plane with the track fitting, skip first 10 mm
    float shr_tkfit_gap10_dedx_V;
    float shr_tkfit_gap10_dedx_U;
    
    unsigned int   shr_tkfit_gap10_nhits_Y; // Number of hits in the 1x4 cm box on the Y plane with the track fitting, skip first 10 mm
    unsigned int   shr_tkfit_gap10_nhits_V;
    unsigned int   shr_tkfit_gap10_nhits_U;
    
    float CylFrac1h_1cm;  // Frac of spacepoints of the leading shower within 1cm of the shower axis. Only in the first half of the shower
    float CylFrac1h_2cm;  // Frac of spacepoints of the leading shower within 2cm of the shower axis. Only in the first half of the shower
    float CylFrac1h_3cm;  // Frac of spacepoints of the leading shower within 3cm of the shower axis. Only in the first half of the shower
    float CylFrac1h_4cm;  // Frac of spacepoints of the leading shower within 4cm of the shower axis. Only in the first half of the shower
    float CylFrac1h_5cm;  // Frac of spacepoints of the leading shower within 5cm of the shower axis. Only in the first half of the shower

    float CylFrac2h_1cm;  // Frac of spacepoints of the leading shower within 1cm of the shower axis. Only in the second half of the shower
    float CylFrac2h_2cm;  // Frac of spacepoints of the leading shower within 2cm of the shower axis. Only in the second half of the shower
    float CylFrac2h_3cm;  // Frac of spacepoints of the leading shower within 3cm of the shower axis. Only in the second half of the shower
    float CylFrac2h_4cm;  // Frac of spacepoints of the leading shower within 4cm of the shower axis. Only in the second half of the shower
    float CylFrac2h_5cm;  // Frac of spacepoints of the leading shower within 5cm of the shower axis. Only in the second half of the shower

    float CylFrac_1cm;  // Frac of spacepoints of the leading shower within 1cm of the shower axis.
    float CylFrac_2cm;  // Frac of spacepoints of the leading shower within 2cm of the shower axis.
    float CylFrac_3cm;  // Frac of spacepoints of the leading shower within 3cm of the shower axis.
    float CylFrac_4cm;  // Frac of spacepoints of the leading shower within 4cm of the shower axis.
    float CylFrac_5cm;  // Frac of spacepoints of the leading shower within 5cm of the shower axis.

    float DeltaMed;   // to set branch
    float DeltaMed1h; // to set branch
    float DeltaMed2h; // to set branch

    float DeltaRMS;   // to set branch
    float DeltaRMS1h; // to set branch
    float DeltaRMS2h; // 0 to 20

    float shrPCA1CMed_5cm; // 0 to 1

    float shrMCSMom; // 0 to 200


    float shr_chipr;            // Shower: Chi2 proton score for the leading shower (with the shower reconstructed as track)
    float shr_chimu;            // Shower: Chi2 muon score for the leading shower (with the shower reconstructed as track)
    float shr_bragg_p;          // Shower: Proton Bragg likelihood score for the leading shower (with the shower reconstructed as track)
    float shr_bragg_mu;         // Shower: Muon Bragg likelihood score for the leading shower (with the shower reconstructed as track)
    float shr_bragg_mip;        // Shower: MIP Bragg likelihood score for the leading shower (with the shower reconstructed as track)
    float shr_bragg_kaon;       // Shower: Kaon Bragg likelihood for the leading shower (with the shower reconstructed as track)
    float shr_bragg_pion;       // Shower: Pion Bragg likelihood for the leading shower (with the shower reconstructed as track)
    
    float tksh_distance;        // Distance between leading shower vertex and longest track vertex
    float tksh_angle;           // Angle between leading shower vertex and longest track vertex

    float shr_distance;         // Shower: Distance between leading shower vertex and reconstructed neutrino vertex. Labelled as shower_vtx_dist in technote
    float shr_score;            // Shower: Pandora track score for the leading shower
    
    int   shr_bkt_pdg;          // Shower: PDG code of the MCParticle matched to the leading shower
    float shr_bkt_purity;       // Shower: Purity of the leading shower
    float shr_bkt_completeness; // Shower: Completeness of the leading shower
    float shr_bkt_E;            // Shower: Energy of the MCParticle matched to the leading shower
    
    float trk_len;              // Track: Length of the longest track
    float trk_theta;            // Track: Reconstructed theta angle for the longest track
    float trk_phi;              // Track: Reconstructed phi angle for the longest track
    float trk_energy;           // Track: Energy of the longest track assuming it's a proton and using stopping power in LAr 
    float trk_energy_muon;      // Track: Energy of the longest track assuming it's a muon and using stopping power in LAr
    float trk_energy_muon_mcs;  // Track: Energy of the longest track assuming it's a muon and using MCS
    float trk_energy_tot;       // Track: Sum of the track energies assuming they are protons and using stopping power in LAr
    float trk_energy_muon_tot;  // Track: Sum of the track energies assuming they are muons and using stopping power in LAr
    float trk_distance;         // Track: Distance between longest track and reconstructed neutrino vertex
    float trk_score;            // Track: Pandora track score for the longest track
    int   trk_bkt_pdg;          // Track: Backtrack
    float trk_bkt_purity;       // Track: Backtrack Purity
    float trk_bkt_completeness; // Track: Backtrack Completeness
    float trk_bkt_E;            // Track: Backtrack Energy
    float trk_chipr_best;       // Track: Chi Prob Best
    float trk_chipr_worst;      // Track: Chi Prob Worst
    float trk_chimu_best;       // Track: Chi Mu Best
    float trk_chimu_worst;      // Track: Chi Mu Worst
    float trk_chipr;            // Track: Chi Probabilty
    float trk_chimu;            // Track: Chi Mu
    float trk_pida;             // Track: PIDA score for the longest track
    
    float trk_bragg_p;          // Track Bragg Likelihood P
    float trk_bragg_mu;         // Track Bragg Likelihood Mu
    float trk_bragg_mip;        // Track Bragg Likelihood MIPP
    float trk_bragg_kaon;       // Track Bragg Likelihood Kaon
    float trk_bragg_pion;       // Track Bragg Likelihood Pion
    
    unsigned int   trk_hits_max;         // Track:  Number of hits of the leading track
    unsigned int   shr_hits_max;         // Shower: Number of hits of the leading shower 
    
    float trkshrhitdist0;       // distance between hits of shower and track in 2D on each palne based on hit-hit distances
    float trkshrhitdist1;      
    float trkshrhitdist2;      
    
    int   total_hits_y;         // Total number of hits on the Y plane
    float extra_energy_y;       // Total energy of the unclustered hits on the Y plane
    float trk_energy_hits_tot;  // Track: Sum of the energy of the tracks obtained with the deposited charge
    
    unsigned int   shrsubclusters0;      // Shower: in how many sub-clusters can the shower be broken based on proximity clustering?
    unsigned int   shrsubclusters1;      // Plane 1
    unsigned int   shrsubclusters2;      // Plane 2
    
    float shrclusfrac0;         // what fraction of the total charge does the dominant shower sub-cluster carry?
    float shrclusfrac1;
    float shrclusfrac2;
    
    float shrclusdir0;          // Shower: 2D charge-weighted direction of shower hits calculated from neutrino vertex w.r.t. vertical in plane
    float shrclusdir1; 
    float shrclusdir2;
    
    unsigned int   shr_hits_tot;         // Shower: Total number of shower hits (all showers)
    unsigned int   shr_hits_y_tot;       // Shower: Total number of shower hits on the Y Plane 2
    unsigned int   shr_hits_u_tot;       // Shower: Total number of shower hits on the U Plane 1 
    unsigned int   shr_hits_v_tot;       // Shower: Total number of shower hits on the V Plane 0
    
    unsigned int   trk_hits_tot;         // Track: Total number of track hits
    unsigned int   trk_hits_y_tot;       // Track: Total number of track hits on the Y Plane 2
    unsigned int   trk_hits_u_tot;       // Track: Total number of track hits on the U Plane 0
    unsigned int   trk_hits_v_tot;       // Track: Total number of track hits on the V Plane 1
    
    float elecclusters_U_charge;
    float elecclusters_V_charge;
    float elecclusters_Y_charge;
    int   elecclusters_U_N;
    int   elecclusters_V_N;
    int   elecclusters_Y_N;
    
    unsigned int   n_tracks_contained;    // Reco Number of tracks fully contained in the fiducial volume
    unsigned int   n_showers_contained;   // Reco Number of showers with a starting point within the fiducial volume.
    float matched_E;             // Total kinetic energy of the MCParticles matched to PFParticles
    float hits_ratio;            // Reco Ratio between hits from showers and total number of hits in the slice.
    float contained_fraction;    // Reco Hits in PFParticles contained in the fiducial volume over the total number of clustered hits in the slice.
    float sps_contained_fraction;
    float pt;                    // Total reconstructed transverse momentum, assuming all the tracks are protons and all the showers are electrons
    float p;                     // Total reconstructed momentum, assuming all the tracks are protons and all the showers are electrons
    float pt_assume_muon;        // Total reconstructed transverse momentum, assuming all the tracks are muons and all the showers are electrons
    float p_assume_muon;         // Total reconstructed momentum, assuming all the tracks are muons and all the showers are electrons
    float dvtx;                  // smallest distance between vertex and any boundary
    float dtrk;                  // smallest distance between any track start/end point and any boundary
    float contained_sps_ratio;
    
    float CosmicIP;              // Reco Closest distance between shower start and space points associated to tracks flagged as cosmics.
    float CosmicIPAll3D;         // 3D distance of shower start from closest spacepoint of any pfp not in the neutrino slice
    float CosmicDirAll3D;        // cosine of 3D direction difference between shower and closest pfp not in the neutrino slice
    float CosmicIPAll2DEnds;     // 2D distance of shower cluster endpoints from closest cluster endpoints from any pfp not in the neutrino slice (minimum fron all planes)
    float CosmicDirAll2DEnds;    // cosine of 2D direction difference between shower cluster and closest cluster from pfps not in the neutrino slice
    float CosmicIPAll2DOvlp;     // 2D distance of shower cluster endpoints from closest line connecting cluster endpoints from any pfp not in the neutrino slice (minimum fron all planes)
    float CosmicDirAll2DOvlp;    // cosine of 2D direction difference between shower cluster and cluster whose endpoints for the closest line 
    
    float leeweight;
    float weightSplineTimesTune;
    float ppfx_cv;    // Weight from PPFX CV
    
    float true_pt;         // Total Pt of all MC particles
    float true_pt_visible; // Total Visible Pt of all MC particles
    float true_p;          // True momentum of all MC particles
    float true_p_visible;  // True visible momentum of all MC particles
    float true_e_visible;  // True visible energy of all MC particles
    
    float opfilter_pe_beam; // Common Optical Filter (beam window)
    float opfilter_pe_veto; // Common Optical Filter (michel veto)
    
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
    float true_nu_px;            // True Neutrino Px
    float true_nu_py;            // True Neutrino Py
    float true_nu_pz;            // True Neutrino Pz
    
    float reco_nu_vtx_x;         // Reco Neutrino Vtx x
    float reco_nu_vtx_y;         // Reco Neutrino Vtx y
    float reco_nu_vtx_z;         // Reco Neutrino Vtx z
    float reco_nu_vtx_sce_x;     // Reco Reconstructed neutrino interaction vertex in (x,y,z) coordinates. The space charged correction is applied
    float reco_nu_vtx_sce_y;     // Reco Neutrino Vtx Space Charge y
    float reco_nu_vtx_sce_z;     // Reco Neutrino Vtx Space Charge z
    
    int   nmuon;   // Truth Number of Muons
    float muon_e;  // Truth Muon Energy
    float muon_c;  // Muon Completeness
    float muon_p;  // Muon Purity
     
    int   nelec;   // Truth Number of electrons
    float elec_e;  // Truth Electron Energy
    float elec_c;  // Electron Completeness
    float elec_p;  // Electron Purity
    
    float elec_vx;   // Truth Electron Vtx x
    float elec_vy;   // Truth Electron Vtx y
    float elec_vz;   // Truth Electron Vtx z
    float elec_px;   // Truth Electron Px
    float elec_py;   // Truth Electron Py
    float elec_pz;   // Truth Electron Pz

    int   npi0;     // Truth Number of Pi0
    float pi0_e;    // Truth Pi0 Energy
    float pi0_c;    // Pi0 Completeness
    float pi0_p;    // Pi0 Purity
    
    int   nneutron; // Truth Number of Neutrons
    
    int   nproton;  // Truth Number of Protons
    float proton_e; // Truth Proton Energy
    float proton_c; // Proton Completeness
    float proton_p; // Proton Purity
    
    int   npion;    // Truth Number of Pions
    float pion_e;   // Truth Pion Energy
    float pion_c;   // Pion Completeness
    float pion_p;   // Pion Purity
    
    int   nslice{0};     // Reco Number of neutrino slices identified by the SliceID. Values are 0 or 1.
    int   crtveto;    // Reco Boolean variable checking if the event passes the CRT veto.
    float crthitpe;   // CRT Hit PE
    int   category;   // Truth Category
    float lep_e;
    int   pass;
    int   swtrig;     // Software Trigger
    int   swtrig_pre;     // Software Trigger before change
    int   swtrig_post;     // Software Trigger after change
    int   evnhits;    // Reco Number of hits in the event
    int   slpdg;      // Reco Slice PDG of primary PFP
    int   slnhits;    // Reco Slice Number of Hits within it
    int   n_pfps;     // Reco Number of pfp's
    int   n_tracks;   // Reco Number of Tracks identified with a track score cut
    int   n_showers;  // Reco Number of Showers with a shower score cut
    
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
    
    float NeutrinoEnergy0;  // Reco Neutrino energy from anab::calorimetry
    float NeutrinoEnergy1;
    float NeutrinoEnergy2;
    
    float SliceCaloEnergy0; // Reco Total slice energy from anab::calorimetry
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
    
    float nu_completeness_from_pfp; // Neutrino Completeness from PFP (how many of the hits reconstructed for the neutrino were from the true neutrino? )
    float nu_purity_from_pfp;       // Neutrino Purity from PFP (how many out of all the hits are the neutrino) (needs to be higher than 50% otherwise tag as a mixed)
    int   n_tracks_pandora;         // Number of Tracks Returned by Pandora
    
    double _closestNuCosmicDist; // Distance between the neutrino vertex and (closest?) cosmic trajectory tagged from CRT

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

    std::vector<float> *pfp_generation_v          = nullptr; // Vec of PFP generation. 1 is primary
    std::vector<float> *pfp_trk_daughters_v       = nullptr; // Vec PFP Track Daughters
    std::vector<float> *pfp_shr_daughters_v       = nullptr; // Vec PFP Shower Daughters
    std::vector<float> *trk_score_v               = nullptr; // Vec PFP track score (must also be > 0)
    std::vector<float> *pfpdg_v                   = nullptr; // Vec PFP in the slice
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
    
    std::vector<int>   *mc_pdg_v  = nullptr;  // True: Vector of all MC particles
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
    
    TH1D* h_2D_CV_UW_PPFX_ratio_nue;
    TH1D* h_2D_CV_UW_PPFX_ratio_nuebar;
    TH1D* h_2D_CV_UW_PPFX_ratio_numu;
    TH1D* h_2D_CV_UW_PPFX_ratio_numubar;

};

#endif
