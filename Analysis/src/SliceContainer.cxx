#include "../include/SliceContainer.h"

// -----------------------------------------------------------------------------
template<typename T> void SliceContainer::Initialise(T *tree, int type, Utility util){

    std::cout << "Initalising Slice Container" << std::endl;
    _util = util;

    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("sub", &sub);
    tree->SetBranchAddress("evt", &evt);
    
    tree->SetBranchAddress("shr_energy_tot", &shr_energy_tot);
    tree->SetBranchAddress("shr_energy", &shr_energy);
    tree->SetBranchAddress("shr_energy_tot_cali", &shr_energy_tot_cali);
    tree->SetBranchAddress("shr_energy_cali", &shr_energy_cali);
    tree->SetBranchAddress("shr_theta", &shr_theta);
    tree->SetBranchAddress("shr_phi", &shr_phi);
    tree->SetBranchAddress("shr_pca_0", &shr_pca_0);
    tree->SetBranchAddress("shr_pca_1", &shr_pca_1);
    tree->SetBranchAddress("shr_pca_2", &shr_pca_2);
    tree->SetBranchAddress("shr_px", &shr_px);
    tree->SetBranchAddress("shr_py", &shr_py);
    tree->SetBranchAddress("shr_pz", &shr_pz);
    tree->SetBranchAddress("shr_openangle", &shr_openangle);
    tree->SetBranchAddress("shr_tkfit_start_x", &shr_tkfit_start_x);
    tree->SetBranchAddress("shr_tkfit_start_y", &shr_tkfit_start_y);
    tree->SetBranchAddress("shr_tkfit_start_z", &shr_tkfit_start_z);
    tree->SetBranchAddress("shr_tkfit_theta", &shr_tkfit_theta);
    tree->SetBranchAddress("shr_tkfit_phi", &shr_tkfit_phi);
    tree->SetBranchAddress("shr_start_x", &shr_start_x);
    tree->SetBranchAddress("shr_start_y", &shr_start_y);
    tree->SetBranchAddress("shr_start_z", &shr_start_z);
    tree->SetBranchAddress("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y);
    tree->SetBranchAddress("shr_tkfit_dedx_V", &shr_tkfit_dedx_V);
    tree->SetBranchAddress("shr_tkfit_dedx_U", &shr_tkfit_dedx_U);
    tree->SetBranchAddress("shr_tkfit_nhits_Y", &shr_tkfit_nhits_Y);
    tree->SetBranchAddress("shr_tkfit_nhits_V", &shr_tkfit_nhits_V);
    tree->SetBranchAddress("shr_tkfit_nhits_U", &shr_tkfit_nhits_U);
    tree->SetBranchAddress("shr_trkfitmedangle", &shr_trkfitmedangle);
    tree->SetBranchAddress("shrmoliereavg", &shrmoliereavg);
    tree->SetBranchAddress("shrmoliererms", &shrmoliererms);
    
    tree->SetBranchAddress("ismerged",       &ismerged);
    tree->SetBranchAddress("merge_bestdot",  &merge_bestdot);
    tree->SetBranchAddress("merge_bestdist", &merge_bestdist);
    tree->SetBranchAddress("merge_vtx_x",    &merge_vtx_x);
    tree->SetBranchAddress("merge_vtx_y",    &merge_vtx_y);
    tree->SetBranchAddress("merge_vtx_z",    &merge_vtx_z);
    // tree->SetBranchAddress("merge_tk_ipfp", &merge_tk_ipfp);
    
    tree->SetBranchAddress("shr_tkfit_2cm_dedx_Y",    &shr_tkfit_2cm_dedx_Y);
    tree->SetBranchAddress("shr_tkfit_2cm_dedx_V",    &shr_tkfit_2cm_dedx_V);
    tree->SetBranchAddress("shr_tkfit_2cm_dedx_U",    &shr_tkfit_2cm_dedx_U);
    tree->SetBranchAddress("shr_tkfit_2cm_nhits_Y",   &shr_tkfit_2cm_nhits_Y  );
    tree->SetBranchAddress("shr_tkfit_2cm_nhits_V",   &shr_tkfit_2cm_nhits_V  );
    tree->SetBranchAddress("shr_tkfit_2cm_nhits_U",   &shr_tkfit_2cm_nhits_U  );
    tree->SetBranchAddress("shr_tkfit_gap05_dedx_Y",  &shr_tkfit_gap05_dedx_Y );
    tree->SetBranchAddress("shr_tkfit_gap05_dedx_V",  &shr_tkfit_gap05_dedx_V );
    tree->SetBranchAddress("shr_tkfit_gap05_dedx_U",  &shr_tkfit_gap05_dedx_U );
    tree->SetBranchAddress("shr_tkfit_gap05_nhits_Y", &shr_tkfit_gap05_nhits_Y );
    tree->SetBranchAddress("shr_tkfit_gap05_nhits_V", &shr_tkfit_gap05_nhits_V );
    tree->SetBranchAddress("shr_tkfit_gap05_nhits_U", &shr_tkfit_gap05_nhits_U );
    tree->SetBranchAddress("shr_tkfit_gap10_dedx_Y",  &shr_tkfit_gap10_dedx_Y );
    tree->SetBranchAddress("shr_tkfit_gap10_dedx_V",  &shr_tkfit_gap10_dedx_V );
    tree->SetBranchAddress("shr_tkfit_gap10_dedx_U",  &shr_tkfit_gap10_dedx_U );
    tree->SetBranchAddress("shr_tkfit_gap10_nhits_Y", &shr_tkfit_gap10_nhits_Y );
    tree->SetBranchAddress("shr_tkfit_gap10_nhits_V", &shr_tkfit_gap10_nhits_V );
    tree->SetBranchAddress("shr_tkfit_gap10_nhits_U", &shr_tkfit_gap10_nhits_U );

    tree->SetBranchAddress("all_shr_hits",     &all_shr_hits);
    tree->SetBranchAddress("all_shr_energies", &all_shr_energies);

    tree->SetBranchAddress("shrPCA1CMed_5cm", &shrPCA1CMed_5cm );

    tree->SetBranchAddress("shrMCSMom", &shrMCSMom );

    tree->SetBranchAddress("shr_chipr",      &shr_chipr);
    tree->SetBranchAddress("shr_chimu",      &shr_chimu);
    tree->SetBranchAddress("shr_bragg_p",    &shr_bragg_p);
    tree->SetBranchAddress("shr_bragg_mu",   &shr_bragg_mu);
    tree->SetBranchAddress("shr_bragg_mip",  &shr_bragg_mip);
    tree->SetBranchAddress("shr_bragg_kaon", &shr_bragg_kaon);
    tree->SetBranchAddress("shr_bragg_pion", &shr_bragg_pion);
    
    tree->SetBranchAddress("tksh_distance", &tksh_distance);
    tree->SetBranchAddress("tksh_angle",    &tksh_angle);
    
    tree->SetBranchAddress("shr_distance", &shr_distance);
    tree->SetBranchAddress("shr_score", &shr_score);
    tree->SetBranchAddress("shr_bkt_pdg", &shr_bkt_pdg);
    tree->SetBranchAddress("shr_bkt_purity", &shr_bkt_purity);
    tree->SetBranchAddress("shr_bkt_completeness", &shr_bkt_completeness);
    tree->SetBranchAddress("shr_bkt_E", &shr_bkt_E);
    
    tree->SetBranchAddress("trk_len", &trk_len);
    tree->SetBranchAddress("trk_theta", &trk_theta);
    tree->SetBranchAddress("trk_phi", &trk_phi);
    tree->SetBranchAddress("trk_energy", &trk_energy);
    tree->SetBranchAddress("trk_energy_muon", &trk_energy_muon);
    tree->SetBranchAddress("trk_energy_muon_mcs", &trk_energy_muon_mcs);
    tree->SetBranchAddress("trk_energy_tot", &trk_energy_tot);
    tree->SetBranchAddress("trk_energy_muon_tot", &trk_energy_muon_tot);
    tree->SetBranchAddress("trk_distance", &trk_distance);
    tree->SetBranchAddress("trk_score", &trk_score);
    tree->SetBranchAddress("trk_bkt_pdg", &trk_bkt_pdg);
    tree->SetBranchAddress("trk_bkt_purity", &trk_bkt_purity);
    tree->SetBranchAddress("trk_bkt_completeness", &trk_bkt_completeness);
    tree->SetBranchAddress("trk_bkt_E", &trk_bkt_E);
    tree->SetBranchAddress("trk_chipr_best", &trk_chipr_best);
    tree->SetBranchAddress("trk_chipr_worst", &trk_chipr_worst);
    tree->SetBranchAddress("trk_chimu_best", &trk_chimu_best);
    tree->SetBranchAddress("trk_chimu_worst", &trk_chimu_worst);
    tree->SetBranchAddress("trk_chipr", &trk_chipr);
    tree->SetBranchAddress("trk_chimu", &trk_chimu);
    tree->SetBranchAddress("trk_pida", &trk_pida);
    tree->SetBranchAddress("trk_bragg_p", &trk_bragg_p);
    tree->SetBranchAddress("trk_bragg_mu", &trk_bragg_mu);
    tree->SetBranchAddress("trk_bragg_mip", &trk_bragg_mip);
    tree->SetBranchAddress("trk_bragg_kaon", &trk_bragg_kaon);
    tree->SetBranchAddress("trk_bragg_pion", &trk_bragg_pion);
    // tree->SetBranchAddress("trk_hits_max", &trk_hits_max);
    tree->SetBranchAddress("shr_hits_max", &shr_hits_max);
    tree->SetBranchAddress("trkshrhitdist0", &trkshrhitdist0);
    tree->SetBranchAddress("trkshrhitdist1", &trkshrhitdist1);
    tree->SetBranchAddress("trkshrhitdist2", &trkshrhitdist2);
    // tree->SetBranchAddress("total_hits_y", &total_hits_y);
    tree->SetBranchAddress("extra_energy_y", &extra_energy_y);
    tree->SetBranchAddress("trk_energy_hits_tot", &trk_energy_hits_tot);
    tree->SetBranchAddress("shrsubclusters0", &shrsubclusters0);
    tree->SetBranchAddress("shrsubclusters1", &shrsubclusters1);
    tree->SetBranchAddress("shrsubclusters2", &shrsubclusters2);
    
    tree->SetBranchAddress("shrclusfrac0", &shrclusfrac0);
    tree->SetBranchAddress("shrclusfrac1", &shrclusfrac1);
    tree->SetBranchAddress("shrclusfrac2", &shrclusfrac2);
    tree->SetBranchAddress("shrclusdir0", &shrclusdir0);
    tree->SetBranchAddress("shrclusdir1", &shrclusdir1);
    tree->SetBranchAddress("shrclusdir2", &shrclusdir2);
    tree->SetBranchAddress("shr_hits_tot", &shr_hits_tot);
    tree->SetBranchAddress("shr_hits_y_tot", &shr_hits_y_tot);
    tree->SetBranchAddress("shr_hits_u_tot", &shr_hits_u_tot);
    tree->SetBranchAddress("shr_hits_v_tot", &shr_hits_v_tot);
    tree->SetBranchAddress("trk_hits_tot", &trk_hits_tot);
    tree->SetBranchAddress("trk_hits_y_tot", &trk_hits_y_tot);
    tree->SetBranchAddress("trk_hits_u_tot", &trk_hits_u_tot);
    tree->SetBranchAddress("trk_hits_v_tot", &trk_hits_v_tot);
    tree->SetBranchAddress("n_tracks_contained",  &n_tracks_contained);
    tree->SetBranchAddress("n_showers_contained", &n_showers_contained);
    
    tree->SetBranchAddress("matched_E", &matched_E);
    tree->SetBranchAddress("hits_ratio", &hits_ratio);
    tree->SetBranchAddress("contained_fraction", &contained_fraction);
    tree->SetBranchAddress("sps_contained_fraction", &sps_contained_fraction);
    tree->SetBranchAddress("pt", &pt);
    tree->SetBranchAddress("p", &p);
    tree->SetBranchAddress("pt_assume_muon", &pt_assume_muon);
    tree->SetBranchAddress("p_assume_muon", &p_assume_muon);
    tree->SetBranchAddress("dvtx", &dvtx);
    tree->SetBranchAddress("dtrk", &dtrk);
    tree->SetBranchAddress("contained_sps_ratio", &contained_sps_ratio);
    tree->SetBranchAddress("CosmicIP", &CosmicIP);
    tree->SetBranchAddress("CosmicIPAll3D", &CosmicIPAll3D);
    tree->SetBranchAddress("CosmicDirAll3D", &CosmicDirAll3D);
    tree->SetBranchAddress("CosmicIPAll2DEnds", &CosmicIPAll2DEnds);
    tree->SetBranchAddress("CosmicDirAll2DEnds", &CosmicDirAll2DEnds);
    tree->SetBranchAddress("CosmicIPAll2DOvlp", &CosmicIPAll2DOvlp);
    tree->SetBranchAddress("CosmicDirAll2DOvlp", &CosmicDirAll2DOvlp);
    tree->SetBranchAddress("leeweight", &leeweight);
    tree->SetBranchAddress("true_pt", &true_pt);
    tree->SetBranchAddress("true_pt_visible", &true_pt_visible);
    tree->SetBranchAddress("true_p", &true_p);
    tree->SetBranchAddress("true_p_visible", &true_p_visible);
    tree->SetBranchAddress("true_e_visible", &true_e_visible);
    
    tree->SetBranchAddress("_opfilter_pe_beam", &opfilter_pe_beam);
    tree->SetBranchAddress("_opfilter_pe_veto", &opfilter_pe_veto);
    
    tree->SetBranchAddress("nu_pdg", &nu_pdg);
    tree->SetBranchAddress("ccnc",   &ccnc);
    tree->SetBranchAddress("interaction", &interaction);
    tree->SetBranchAddress("nu_e",   &nu_e);
    tree->SetBranchAddress("nu_pt",  &nu_pt);
    tree->SetBranchAddress("theta",  &theta);
    tree->SetBranchAddress("isVtxInFiducial",  &isVtxInFiducial);
    tree->SetBranchAddress("truthFiducial",    &truthFiducial);
    tree->SetBranchAddress("reco_e",           &reco_e);
    
    tree->SetBranchAddress("true_nu_vtx_t", &true_nu_vtx_t);
    tree->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x);
    tree->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y);
    tree->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z);
    tree->SetBranchAddress("true_nu_vtx_sce_x", &true_nu_vtx_sce_x);
    tree->SetBranchAddress("true_nu_vtx_sce_y", &true_nu_vtx_sce_y);
    tree->SetBranchAddress("true_nu_vtx_sce_z", &true_nu_vtx_sce_z);
    tree->SetBranchAddress("true_nu_px", &true_nu_px);
    tree->SetBranchAddress("true_nu_py", &true_nu_py);
    tree->SetBranchAddress("true_nu_pz", &true_nu_pz);
    
    tree->SetBranchAddress("reco_nu_vtx_x", &reco_nu_vtx_x);
    tree->SetBranchAddress("reco_nu_vtx_y", &reco_nu_vtx_y);
    tree->SetBranchAddress("reco_nu_vtx_z", &reco_nu_vtx_z);
    tree->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x);
    tree->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y);
    tree->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z);
    
    tree->SetBranchAddress("nmuon", &nmuon);
    tree->SetBranchAddress("muon_e", &muon_e);
    tree->SetBranchAddress("muon_c", &muon_c);
    tree->SetBranchAddress("muon_p", &muon_p);
    tree->SetBranchAddress("nelec", &nelec);
    
    tree->SetBranchAddress("elec_e", &elec_e);
    tree->SetBranchAddress("elec_c", &elec_c);
    tree->SetBranchAddress("elec_p", &elec_p);
    tree->SetBranchAddress("elec_vx", &elec_vx);
    tree->SetBranchAddress("elec_vy", &elec_vy);
    tree->SetBranchAddress("elec_vz", &elec_vz);
    tree->SetBranchAddress("elec_px", &elec_px);
    tree->SetBranchAddress("elec_py", &elec_py);
    tree->SetBranchAddress("elec_pz", &elec_pz);
    
    tree->SetBranchAddress("npi0", &npi0);
    tree->SetBranchAddress("pi0_e", &pi0_e);
    tree->SetBranchAddress("pi0_c", &pi0_c);
    tree->SetBranchAddress("pi0_p", &pi0_p);
    tree->SetBranchAddress("nneutron", &nneutron);
    tree->SetBranchAddress("nproton", &nproton);
    
    tree->SetBranchAddress("proton_e", &proton_e);
    tree->SetBranchAddress("proton_c", &proton_c);
    tree->SetBranchAddress("proton_p", &proton_p);
    
    tree->SetBranchAddress("npion", &npion);
    tree->SetBranchAddress("pion_e", &pion_e);
    tree->SetBranchAddress("pion_c", &pion_c);
    tree->SetBranchAddress("pion_p", &pion_p);
    tree->SetBranchAddress("nslice", &nslice);
    tree->SetBranchAddress("crtveto", &crtveto);
    tree->SetBranchAddress("crthitpe", &crthitpe);
    tree->SetBranchAddress("category", &category);
    tree->SetBranchAddress("lep_e", &lep_e);
    tree->SetBranchAddress("pass", &pass);
    
    tree->SetBranchAddress("swtrig", &swtrig);
    tree->SetBranchAddress("swtrig_pre", &swtrig_pre);
    tree->SetBranchAddress("swtrig_post", &swtrig_post);
    
    tree->SetBranchAddress("evnhits", &evnhits);
    tree->SetBranchAddress("slpdg", &slpdg);
    tree->SetBranchAddress("slnhits", &slnhits);
    tree->SetBranchAddress("n_pfps", &n_pfps);
    tree->SetBranchAddress("n_tracks", &n_tracks);
    tree->SetBranchAddress("n_showers", &n_showers);

    tree->SetBranchAddress("topological_score", &topological_score);
    tree->SetBranchAddress("slclustfrac", &slclustfrac);
    tree->SetBranchAddress("endmuonmichel", &endmuonmichel);
    tree->SetBranchAddress("filter_antibdt", &filter_antibdt);
    tree->SetBranchAddress("filter_ncpi0", &filter_ncpi0);
    tree->SetBranchAddress("filter_pi0", &filter_pi0);
    tree->SetBranchAddress("filter_ccinclusive", &filter_ccinclusive);
    tree->SetBranchAddress("flash_pe", &flash_pe);
    tree->SetBranchAddress("flash_time", &flash_time);
    tree->SetBranchAddress("nu_flashmatch_score", &nu_flashmatch_score);
    tree->SetBranchAddress("best_cosmic_flashmatch_score", &best_cosmic_flashmatch_score);
    
    tree->SetBranchAddress("NeutrinoEnergy0", &NeutrinoEnergy0);
    tree->SetBranchAddress("NeutrinoEnergy1", &NeutrinoEnergy1);
    tree->SetBranchAddress("NeutrinoEnergy2", &NeutrinoEnergy2);
    
    tree->SetBranchAddress("SliceCaloEnergy0", &SliceCaloEnergy0);
    tree->SetBranchAddress("SliceCaloEnergy1", &SliceCaloEnergy1);
    tree->SetBranchAddress("SliceCaloEnergy2", &SliceCaloEnergy2);
    
    tree->SetBranchAddress("nnoise_pl1", &nnoise_pl1);
    tree->SetBranchAddress("nslhits_pl1", &nslhits_pl1);
    tree->SetBranchAddress("nslnoise_pl1", &nslnoise_pl1);
    tree->SetBranchAddress("nhits_pl1", &nhits_pl1);
    tree->SetBranchAddress("frac_slnoise_pl1", &frac_slnoise_pl1);

    tree->SetBranchAddress("evnunhits", &evnunhits);
    tree->SetBranchAddress("evlepnhits", &evlepnhits);
    tree->SetBranchAddress("evpronhits", &evpronhits);
    tree->SetBranchAddress("evpi1nhits", &evpi1nhits);
    tree->SetBranchAddress("evpi0nhits", &evpi0nhits);
    tree->SetBranchAddress("evneunhits", &evneunhits);
    tree->SetBranchAddress("evgamnhits", &evgamnhits);
    tree->SetBranchAddress("evothnhits", &evothnhits);
    
    tree->SetBranchAddress("slnunhits", &slnunhits);
    tree->SetBranchAddress("sllepnhits", &sllepnhits);
    tree->SetBranchAddress("slpronhits", &slpronhits);
    tree->SetBranchAddress("slpi1nhits", &slpi1nhits);
    tree->SetBranchAddress("slpi0nhits", &slpi0nhits);
    tree->SetBranchAddress("slneunhits", &slneunhits);
    tree->SetBranchAddress("slgamnhits", &slgamnhits);
    tree->SetBranchAddress("slothnhits", &slothnhits);
    
    tree->SetBranchAddress("nu_completeness_from_pfp", &nu_completeness_from_pfp);
    tree->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp);
    // tree->SetBranchAddress("n_tracks_pandora", &n_tracks_pandora);
    
    if (std::string(_util.run_period) != "1") tree->SetBranchAddress("_closestNuCosmicDist",&_closestNuCosmicDist);
    
    tree->SetBranchAddress("pfp_generation_v",        &pfp_generation_v);
    tree->SetBranchAddress("pfp_trk_daughters_v",     &pfp_trk_daughters_v);
    tree->SetBranchAddress("pfp_shr_daughters_v",     &pfp_shr_daughters_v);
    tree->SetBranchAddress("trk_score_v",             &trk_score_v);
    tree->SetBranchAddress("pfpdg",                   &pfpdg_v);
    tree->SetBranchAddress("pfnhits",                 &pfnhits_v);
    tree->SetBranchAddress("pfnplanehits_U",          &pfnplanehits_U_v);
    tree->SetBranchAddress("pfnplanehits_V",          &pfnplanehits_V_v);
    tree->SetBranchAddress("pfnplanehits_Y",          &pfnplanehits_Y_v);
    tree->SetBranchAddress("pfpplanesubclusters_U",   &pfpplanesubclusters_U_v);
    tree->SetBranchAddress("pfpplanesubclusters_V",   &pfpplanesubclusters_V_v);
    tree->SetBranchAddress("pfpplanesubclusters_Y",   &pfpplanesubclusters_Y_v);
    tree->SetBranchAddress("pfpplanesubhitfracmax_U", &pfpplanesubhitfracmax_U_v);
    tree->SetBranchAddress("pfpplanesubhitfracmax_V", &pfpplanesubhitfracmax_V_v);
    tree->SetBranchAddress("pfpplanesubhitfracmax_Y", &pfpplanesubhitfracmax_Y_v);
    
    tree->SetBranchAddress("mc_pdg",          &mc_pdg_v);
    tree->SetBranchAddress("mc_E",            &mc_E_v);
    tree->SetBranchAddress("mc_vx",           &mc_vx_v);
    tree->SetBranchAddress("mc_vy",           &mc_vy_v);
    tree->SetBranchAddress("mc_vz",           &mc_vz_v);
    tree->SetBranchAddress("mc_endx",         &mc_endx_v);
    tree->SetBranchAddress("mc_endy",         &mc_endy_v);
    tree->SetBranchAddress("mc_endz",         &mc_endz_v);
    tree->SetBranchAddress("mc_px",           &mc_px_v);
    tree->SetBranchAddress("mc_py",           &mc_py_v);
    tree->SetBranchAddress("mc_pz",           &mc_pz_v);
    tree->SetBranchAddress("mc_completeness", &mc_completeness_v);
    tree->SetBranchAddress("mc_purity",       &mc_purity_v);
    
    // MC specific branches
    if (type == _util.k_mc || type == _util.k_dirt){

        tree->SetBranchAddress("weightSplineTimesTune",      &weightSplineTimesTune);
        tree->SetBranchAddress("ppfx_cv",                    &ppfx_cv);
        
        // weightstree->SetBranchAddress("weights", &_mapWeight);
        // tree->SetBranchAddress("weights",         &weights_v);
        
        if (type == _util.k_mc){
            tree->SetBranchAddress("weightsGenie",          &weightsGenie);
            tree->SetBranchAddress("weightsReint",          &weightsReint);
            tree->SetBranchAddress("weightsPPFX",           &weightsPPFX);
            tree->SetBranchAddress("knobRPAup",             &knobRPAup);
            tree->SetBranchAddress("knobRPAdn",             &knobRPAdn);
            tree->SetBranchAddress("knobCCMECup",           &knobCCMECup);
            tree->SetBranchAddress("knobCCMECdn",           &knobCCMECdn);
            tree->SetBranchAddress("knobAxFFCCQEup",        &knobAxFFCCQEup);
            tree->SetBranchAddress("knobAxFFCCQEdn",        &knobAxFFCCQEdn);
            tree->SetBranchAddress("knobVecFFCCQEup",       &knobVecFFCCQEup);
            tree->SetBranchAddress("knobVecFFCCQEdn",       &knobVecFFCCQEdn);
            tree->SetBranchAddress("knobDecayAngMECup",     &knobDecayAngMECup);
            tree->SetBranchAddress("knobDecayAngMECdn",     &knobDecayAngMECdn);
            tree->SetBranchAddress("knobThetaDelta2Npiup",  &knobThetaDelta2Npiup);
            tree->SetBranchAddress("knobThetaDelta2Npidn",  &knobThetaDelta2Npidn);
            tree->SetBranchAddress("knobThetaDelta2NRadup", &knobThetaDelta2NRadup);
            tree->SetBranchAddress("knobThetaDelta2NRaddn", &knobThetaDelta2NRaddn);
            tree->SetBranchAddress("knobRPA_CCQE_Reducedup",&knobRPA_CCQE_Reducedup);
            tree->SetBranchAddress("knobRPA_CCQE_Reduceddn",&knobRPA_CCQE_Reduceddn);
            tree->SetBranchAddress("knobNormCCCOHup",       &knobNormCCCOHup);
            tree->SetBranchAddress("knobNormCCCOHdn",       &knobNormCCCOHdn);
            tree->SetBranchAddress("knobNormNCCOHup",       &knobNormNCCOHup);
            tree->SetBranchAddress("knobNormNCCOHdn",       &knobNormNCCOHdn);
        }
    }
    
    
    tree->SetBranchAddress("cosmic_flashmatch_score_v",&cosmic_flashmatch_score_v);

    tree->SetBranchAddress("pi0_shrscore1",  &pi0_shrscore1);
    tree->SetBranchAddress("pi0_shrscore2",  &pi0_shrscore2);
    tree->SetBranchAddress("pi0_dot1",       &pi0_dot1);
    tree->SetBranchAddress("pi0_dot2",       &pi0_dot2);
    tree->SetBranchAddress("pi0_radlen1",    &pi0_radlen1);
    tree->SetBranchAddress("pi0_radlen2",    &pi0_radlen2);
    tree->SetBranchAddress("pi0_gammadot",   &pi0_gammadot);
    tree->SetBranchAddress("pi0_energy1_Y",  &pi0_energy1_Y);
    tree->SetBranchAddress("pi0_energy2_Y",  &pi0_energy2_Y);
    tree->SetBranchAddress("pi0_dedx1_fit_Y",&pi0_dedx1_fit_Y);
    tree->SetBranchAddress("pi0_mass_Y",     &pi0_mass_Y);
    
    tree->SetBranchAddress("trk_pfp_id_v",        &trk_pfp_id_v);
    tree->SetBranchAddress("trk_sce_start_x_v",        &trk_sce_start_x_v);
    tree->SetBranchAddress("trk_sce_start_y_v",        &trk_sce_start_y_v);
    tree->SetBranchAddress("trk_sce_start_z_v",        &trk_sce_start_z_v);
    tree->SetBranchAddress("trk_sce_end_x_v",          &trk_sce_end_x_v);
    tree->SetBranchAddress("trk_sce_end_y_v",          &trk_sce_end_y_v);
    tree->SetBranchAddress("trk_sce_end_z_v",          &trk_sce_end_z_v);
    tree->SetBranchAddress("trk_distance_v",           &trk_distance_v);
    tree->SetBranchAddress("trk_len_v",                &trk_len_v);
    tree->SetBranchAddress("trk_mcs_muon_mom_v",       &trk_mcs_muon_mom_v);
    tree->SetBranchAddress("trk_range_muon_mom_v",     &trk_range_muon_mom_v);
    tree->SetBranchAddress("trk_llr_pid_score_v",      &trk_llr_pid_score_v);

}
// Explicitly initialise the templates for this function which include TTree and TChain
template void SliceContainer::Initialise<TTree>(TTree *tree, int type, Utility util);
template void SliceContainer::Initialise<TChain>(TChain *tree, int type, Utility util );
// -----------------------------------------------------------------------------
void SliceContainer::SliceClassifier(int type){
    
    // MC Specific classsifications
    if (type == _util.k_mc){

        // If its ske data, then all events should have the data classification
        if (_util.isfakedata){
            classification = std::make_pair("data",_util.k_leg_data);
            return;
        }


        
        bool is_in_fv = _util.in_fv(true_nu_vtx_sce_x, true_nu_vtx_sce_y, true_nu_vtx_sce_z);

        // Out of Fiducial Volume Event
        if (!is_in_fv) {
            
            // Cosmic events including mixed
            if (nu_purity_from_pfp >= 0 && nu_purity_from_pfp < 0.5){
                
                classification = std::make_pair("cosmic",_util.k_cosmic);
                return;
            }
            // Not recontructed
            else if (nu_purity_from_pfp < 0){
                classification = std::make_pair("unmatched",_util.k_unmatched);
                return;
            }
            // Out of Fv events
            else {
                classification = std::make_pair("nu_out_fv",_util.k_nu_out_fv);
                return;
            }
        }
        // In FV event
        else {

            // Charged Current 
            if (ccnc == _util.k_CC){

                // NuMu CC
                if (nu_pdg == 14 || nu_pdg == -14){

                    // Purity is low so return cosmic
                    if (nu_purity_from_pfp <= 0.5) {
                        classification = std::make_pair("cosmic",_util.k_cosmic);
                        return;
                    }
                    
                    if (npi0 > 0) {
                        classification = std::make_pair("numu_cc_pi0", _util.k_numu_cc_pi0); // has a pi0
                        return;
                    }
                    else {
                        classification = std::make_pair("numu_cc",_util.k_numu_cc);
                        return;
                    }

                }
                // Nue CC
                else if (nu_pdg == 12){
                    
                    // purity > 0.5% so signal
                    if (nu_purity_from_pfp > 0.5){
                        classification = std::make_pair("nue_cc",       _util.k_nue_cc);    
                        return;
                    }
                    // Most of the slice from pandora was unmatched hits -- probably a cosmic
                    else if (nu_purity_from_pfp >= 0 && nu_purity_from_pfp <= 0.5){
                        classification = std::make_pair("cosmic_nue",   _util.k_cosmic_nue);
                        return;
                    }
                    // These events were not picked up by pandora at all
                    else {
                        classification = std::make_pair("unmatched_nue",_util.k_unmatched_nue); 
                        return;
                    }

                }
                else if (nu_pdg == -12){
                    
                    // purity > 0.5% so signal
                    if (nu_purity_from_pfp > 0.5) {
                        classification = std::make_pair("nuebar_cc",       _util.k_nuebar_cc); 
                        return;

                    }
                    // Most of the slice from pandora was unmatched hits -- probably a cosmic                               
                    else if (nu_purity_from_pfp >= 0 && nu_purity_from_pfp <= 0.5) {
                        classification = std::make_pair("cosmic_nuebar",   _util.k_cosmic_nuebar);    
                        return;
                    } 
                    // These events were not picked up by pandora at all
                    else {
                        classification = std::make_pair("unmatched_nuebar",_util.k_unmatched_nuebar); 
                        return;
                    }

                }
                // Unknown Neutrino Type
                else {
                    std::cout << "Unknown Neutrino Type..., This will also mess up the efficecy if this occurs!" << std::endl;
                    classification = std::make_pair("unmatched",_util.k_unmatched);
                    return;
                }

            }
            // Neutral Current
            else {

                // Purity is low so return cosmic
                if (nu_purity_from_pfp <= 0.5){
                    classification = std::make_pair("cosmic",_util.k_cosmic);
                    return;
                }

                if (npi0 > 0) {
                    classification = std::make_pair("nc_pi0",_util.k_nc_pi0);
                    return;
                }
                else {
                    classification = std::make_pair("nc",_util.k_nc);
                    return;
                }
            }
        
        } // End if in FV

    }
    // Data
    else if (type == _util.k_data){
        classification = std::make_pair("data",_util.k_leg_data);
        return;
    }
    // EXT
    else if (type == _util.k_ext){
        classification = std::make_pair("ext",_util.k_leg_ext);
        return;
        
    }
    // Dirt
    else if (type == _util.k_dirt){
        classification = std::make_pair("dirt",_util.k_leg_dirt);
        return;
    }
    // What is this type?
    else {
        std::cout << "Got a case we are calling unmatched, this is going to mess up the efficiency in the current way!" << std::endl;
        classification = std::make_pair("unmatched",_util.k_unmatched);
        return;
    }
    
}
// -----------------------------------------------------------------------------
std::string SliceContainer::SliceCategory(){

    if (category == _util.k_pandora_nu_e_other) {
        return "nue_other";

    }
    else if (category == _util.k_pandora_nu_e_cc0pi0p ) {
        return "nu_e_cc0pi0p";

    }
    else if (category == _util.k_pandora_nu_e_cc0pinp ) {
        return "nu_e_cc0pinp";

    }
    else if (category == _util.k_pandora_nu_mu_other) {
        return "nu_mu_other";

    }
    else if (category == _util.k_pandora_nu_mu_pi0 ) {
        return "nu_mu_pi0";

    }
    else if (category == _util.k_pandora_nc) {
        return "nc";

    }
    else if (category == _util.k_pandora_nc_pi0 ) {
        return "nc_pi0";

    }
    else if (category == _util.k_pandora_cosmic) {
        return "cosmic";

    }
    else if (category == _util.k_pandora_outfv) {
        return "outfv";

    }
    else if (category == _util.k_pandora_other) {
        return "other";

    }
    else if (category == _util.k_pandora_data) {
        return "data";
    }
    else {
        std::cout << "Unknown Category type"<< std::endl;
        return "unknown";
    }
}
// -----------------------------------------------------------------------------
void SliceContainer::ParticleClassifier(int type){
    
    // MC Specific classsifications
    if (type == _util.k_mc){

        // Electron
        if (shr_bkt_pdg == 11 || shr_bkt_pdg == -11){
            particle_type = std::make_pair("e",_util.k_electron);
            return;
        }
        // Muon
        else if (shr_bkt_pdg == 13 || shr_bkt_pdg == -13){
            particle_type = std::make_pair("muon",_util.k_muon);
            return;
        }
        // Pion
        else if (shr_bkt_pdg == 211 || shr_bkt_pdg == -211){
            particle_type = std::make_pair("e",_util.k_pion);
            return;
        }
        // Photon 
        else if (shr_bkt_pdg == 22 ){
            particle_type = std::make_pair("photon",_util.k_photon);
            return;
        }
        // Proton
        else if (shr_bkt_pdg == 2212){
            particle_type = std::make_pair("p",_util.k_proton);
            return;
        }
        // Neutron
        else if (shr_bkt_pdg == 2112){
            particle_type = std::make_pair("n",_util.k_neutron);
            return;
        }
        // Kaon
        else if (shr_bkt_pdg == 321 || shr_bkt_pdg == -321 ){
            particle_type = std::make_pair("K",_util.k_kaon);
            return;
        }
        // Other stuff is assumed cosmic
        else {
            particle_type = std::make_pair("cosmic",_util.k_part_cosmic);
            return;
        }


    }
    // Data
    else if (type == _util.k_data){
        particle_type = std::make_pair("data",_util.k_part_data);
        return;
    }
    // EXT
    else if (type == _util.k_ext){
        particle_type = std::make_pair("ext",_util.k_part_ext);
        return;
        
    }
    // Dirt
    else if (type == _util.k_dirt){
        particle_type = std::make_pair("dirt",_util.k_part_dirt);
        return;
    }
    // What is this type?
    else {
        particle_type = std::make_pair("unmatched",_util.k_part_unmatched);
        return;
    }
    
}
// -----------------------------------------------------------------------------
void SliceContainer::SliceInteractionType(int type){

    // Only do this for mc, otherwise return data type
    if (type == _util.k_mc || type == _util.k_dirt){
        std::string nu = "temp";
        std::string CCNC = "temp";

        // Get the nu flavour
        if (nu_pdg == 14 ){
            nu = "numu_";
        }
        else if (nu_pdg == -14){
            nu = "numu_bar_";
        }
        else if (nu_pdg == 12){
            nu = "nue_";
        }
        else if (nu_pdg == -12){
            nu = "nue_bar_";
        }
        else {
            nu = "unknown_";
        }

        // The interaction type
        if (ccnc == _util.k_CC){
            CCNC = "cc_";
        }
        else CCNC = "nc_";


        if (interaction == _util.k_qe) {
            genie_interaction = nu + CCNC + "qe";
            return;

        }
        else if (interaction == _util.k_res ) {
            genie_interaction = nu + CCNC + "res";
            return;

        }
        else if (interaction == _util.k_dis ) {
            genie_interaction = nu + CCNC + "dis";
            return;

        }
        else if (interaction == _util.k_coh) {
            genie_interaction = nu + CCNC + "coh";
            return;

        }
        else if (interaction == _util.k_mec) {
            genie_interaction = nu + CCNC + "mec";
            return;

        }
        else {
            genie_interaction = nu + CCNC + "unknown";
            return;
        }
    }
    else{
        genie_interaction = "data";
        return;
    }



}
// -----------------------------------------------------------------------------
void SliceContainer::Pi0Classifier(int type){

    // Only do this for mc, otherwise return data type
    if (type == _util.k_mc){
        std::string nu = "temp";
        std::string CCNC = "temp";

        // Get the nu flavour
        if (nu_pdg == 14 ||  nu_pdg == -14){
            nu = "numu_";
        }
        else if (nu_pdg == 12){
            nu = "nue_";
        }
        else if (nu_pdg == -12){
            nu = "nue_bar_";
        }
        else {
            nu = "unknown_";
        }

        // The interaction type
        if (ccnc == _util.k_CC){
            CCNC = "cc";
        }
        else CCNC = "nc";
        

        if (npi0 > 0){
            pi0_classification = nu + CCNC + "_pi0";
            return;

        }
        else if (npi0 == 0 ) {
            pi0_classification = nu + CCNC + "";
            return;

        }
        else {
            pi0_classification = nu + CCNC + "unknown";
            return;
        }
    }
    else {
        pi0_classification = "data";
        return;
    }



}
// -----------------------------------------------------------------------------
double SliceContainer::GetPPFXCVWeight(){
    
    double weight = 1.0;

    double nu_theta = _util.GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "beam");

    double xbin{1.0},ybin{1.0};

    if (nu_pdg == 14) {
        xbin = h_2D_CV_UW_PPFX_ratio_numu->GetXaxis()->FindBin(nu_e);
        ybin = h_2D_CV_UW_PPFX_ratio_numu->GetYaxis()->FindBin(nu_theta);
        weight = h_2D_CV_UW_PPFX_ratio_numu->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == -14) {
        xbin = h_2D_CV_UW_PPFX_ratio_numubar->GetXaxis()->FindBin(nu_e);
        ybin = h_2D_CV_UW_PPFX_ratio_numubar->GetYaxis()->FindBin(nu_theta);
        weight = h_2D_CV_UW_PPFX_ratio_numubar->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == 12) {
        xbin = h_2D_CV_UW_PPFX_ratio_nue->GetXaxis()->FindBin(nu_e);
        ybin = h_2D_CV_UW_PPFX_ratio_nue->GetYaxis()->FindBin(nu_theta);
        weight = h_2D_CV_UW_PPFX_ratio_nue->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == -12) {
        xbin = h_2D_CV_UW_PPFX_ratio_nuebar->GetXaxis()->FindBin(nu_e);
        ybin = h_2D_CV_UW_PPFX_ratio_nuebar->GetYaxis()->FindBin(nu_theta);
        weight = h_2D_CV_UW_PPFX_ratio_nuebar->GetBinContent(xbin, ybin);
    }

    // Add some catches to remove unphysical weights
    if (std::isinf(weight))      weight = 1.0; 
    if (std::isnan(weight) == 1) weight = 1.0;
    if (weight > 100)            weight = 1.0;

    // std::cout << nu_theta << "  " << nu_e <<  "  " << weight << std::endl;

    return weight;
}
// -----------------------------------------------------------------------------
double SliceContainer::GetdEdxMax(){

    double dedx_max = -1;

    // We want to also use the dedx when it is defined properly. Sometimes, the plane can have hits but an undefined dedx
    // use the dedx where we get the max number of hits and the dedx > 0
    // int temp_shr_hits_u_tot = shr_hits_u_tot;
    // int temp_shr_hits_v_tot = shr_hits_v_tot;
    // int temp_shr_hits_y_tot = shr_hits_y_tot;

    int temp_shr_hits_u_tot = shr_tkfit_nhits_U; // These variables give a bigger difference in run 1 and run 3
    int temp_shr_hits_v_tot = shr_tkfit_nhits_V;
    int temp_shr_hits_y_tot = shr_tkfit_nhits_Y;

    // If the dedx is undefined, set the hits to zero
    if (shr_tkfit_dedx_U <= 0) temp_shr_hits_u_tot = 0;
    if (shr_tkfit_dedx_V <= 0) temp_shr_hits_v_tot = 0;
    if (shr_tkfit_dedx_Y <= 0) temp_shr_hits_y_tot = 0;


    // Collection plane is the largest
    if (temp_shr_hits_y_tot > temp_shr_hits_u_tot && temp_shr_hits_y_tot > temp_shr_hits_v_tot ){
        dedx_max = shr_tkfit_dedx_Y;
    }
    // V Plane is the largest
    else if (temp_shr_hits_v_tot > temp_shr_hits_u_tot && temp_shr_hits_v_tot > temp_shr_hits_y_tot) {
        dedx_max = shr_tkfit_dedx_V;
        
    }
    // U Plane is the largest
    else if (temp_shr_hits_u_tot > temp_shr_hits_v_tot && temp_shr_hits_u_tot > temp_shr_hits_y_tot){
        dedx_max = shr_tkfit_dedx_U;
        // std::cout << shr_tkfit_dedx_U << " " << shr_tkfit_dedx_V << " " << shr_tkfit_dedx_Y<< "  " << shr_theta*(180/3.14)<< "  " << npi0 <<  std::endl;
        // std::cout << shr_hits_u_tot << " " << shr_hits_v_tot << " " << shr_hits_y_tot<<"\n" <<std::endl;
    }
    // Ok one plane was equal, so need to prioritise planes in preference of y, v, u
    else {

        // If y == any other plane, then y wins
        if (temp_shr_hits_y_tot == temp_shr_hits_u_tot || temp_shr_hits_y_tot == temp_shr_hits_v_tot ){
            dedx_max = shr_tkfit_dedx_Y;
           
        }
        // U == V, ALL Y cases have been used up, so default to v
        else if (temp_shr_hits_u_tot == temp_shr_hits_v_tot ){
            dedx_max = shr_tkfit_dedx_V;
            
        }
        else {
            dedx_max = shr_tkfit_dedx_U;
            
            
        }
    }

    if (dedx_max == -1) {
        std::cout << shr_tkfit_dedx_U << " " << shr_tkfit_dedx_V << " " << shr_tkfit_dedx_Y<< std::endl;
        std::cout << "edge case of dedx comparisons, your logic is flawed!" << std::endl;
    }

    return dedx_max;
    
}
// -----------------------------------------------------------------------------
void SliceContainer::SetCVWeight(double weight){
    cv_weight = weight;
}
// -----------------------------------------------------------------------------
void SliceContainer::SetSignal(){

    if (classification.second == _util.k_nue_cc           || classification.second == _util.k_nuebar_cc ||
        classification.second == _util.k_unmatched_nue    || classification.second == _util.k_cosmic_nue ||
        classification.second == _util.k_unmatched_nuebar || classification.second == _util.k_cosmic_nuebar){
            is_signal = true;
    }
    else 
        is_signal = false;

}
// -----------------------------------------------------------------------------
void SliceContainer::SetTrueElectronThetaPhi(){

    TVector3 vec(elec_px, elec_py, elec_pz);

    elec_theta = vec.Theta() * 180.0/3.14159;
    elec_phi   = vec.Phi()   * 180.0/3.14159;

}
// -----------------------------------------------------------------------------
void SliceContainer::SetNuMIAngularVariables(){


    // --- Effective Angle -- //
    // Calculate the angle between the shower direction and the vector from the target to the nu vtx
    TVector3 shr_dir(shr_px, shr_py, shr_pz); // Shower direction
    shr_dir.Unit();
    
    TVector3 v_targ_uboone(-31387.58422, -3316.402543, -60100.2414);
    TVector3 v_nu_vtx(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z);
    TVector3 v_targ_to_vtx = (-1*v_targ_uboone + v_nu_vtx).Unit(); // -1 because the vector points from uboone to tgt, we need the other way around

    // Set the values
    effective_angle = shr_dir.Angle(v_targ_to_vtx) * 180 / 3.14159;
    cos_effective_angle = std::cos(shr_dir.Angle(v_targ_to_vtx));
    
    // Calculate the effective angle using true variables
    TVector3 nu_dir(true_nu_px, true_nu_py, true_nu_pz); 
    nu_dir.Unit();

    TVector3 v_nu_vtx_true(true_nu_vtx_sce_x, true_nu_vtx_sce_y, true_nu_vtx_sce_z);
    TVector3 v_targ_to_vtx_true = (-1*v_targ_uboone + v_nu_vtx_true).Unit(); // -1 because the vector points from uboone to tgt, we need the other way around

    TVector3 elec_dir(elec_px, elec_py, elec_pz); 
    elec_dir.Unit();
    // true_effective_angle = elec_dir.Angle(v_targ_to_vtx_true) * 180 / 3.14159;
    true_effective_angle = elec_dir.Angle(nu_dir) * 180 / 3.14159; // Use dot product of elec dir to nu dir for this as the true value
    
    // --

    // Momentum of neutrino
    nu_p = std::sqrt(true_nu_px*true_nu_px + true_nu_py*true_nu_py + true_nu_pz*true_nu_pz);

    // Reconstructed Shower momentum
    shr_p = std::sqrt(shr_px*shr_px + shr_py*shr_py + shr_pz*shr_pz);
    
    // Momentum of electron
    elec_mom = std::sqrt(elec_px*elec_px + elec_py*elec_py + elec_pz*elec_pz);

    // --

    TVector3 nu_p_vec(true_nu_px, true_nu_py, true_nu_pz);

    // True nue theta in BNB coordinates (up from beam dir)
    nu_theta = nu_p_vec.Theta() * 180.0/3.14159;
    
    // True nue phi in BNB coordinates (around beam dir)
    nu_phi = nu_p_vec.Phi()   * 180.0/3.14159;
    
    // True nue angle from numi beamline 
    nu_angle = _util.GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "beam"); 

    // True nue angle wrt numi target to uboone vector
    nu_angle_targ = _util.GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "target"); 

    // True electron angle wrt numi target to uboone vector
    elec_ang_targ = _util.GetNuMIAngle(elec_px, elec_py, elec_pz, "target");

    // --- Calculate the dot-product of the proxy for neutrino direction --- //
    // (the vector from the target to the reco nu vtx) and the true nu angle

    reco_true_nu_ang = nu_dir.Angle(v_targ_to_vtx) * 180/3.14159;
    // ---

    // The angle of the reconstructed shower relative to the NuMI target to detector direction
    shr_ang_numi = _util.GetNuMIAngle(shr_px, shr_py, shr_pz, "target");

}
// -----------------------------------------------------------------------------
void SliceContainer::CalibrateShowerEnergy(){

    // Divide the shower energy by 0.83 to calibrate it properly. 
    shr_energy_cali= shr_energy_cali/0.83;

}
// -----------------------------------------------------------------------------
void SliceContainer::SetThresholdEvent(){

    // Below threshold of nue or electron energy
    if (nu_e < _util.energy_threshold || elec_e < _util.elec_threshold){
        if (nu_pdg == 12){
            classification = std::make_pair("thr_nue",_util.k_thr_nue);
        }
        if (nu_pdg == -12){
            classification = std::make_pair("thr_nuebar",_util.k_thr_nuebar);
        }
    }
}