#include "../include/TreeHelper.h"

// -----------------------------------------------------------------------------
void TreeHelper::Initialise(int type, const char* file_out, Utility _utility ){

    std::cout << "Initalising Tree Helper..." << std::endl;

    _util = _utility;

    std::string file_out_str = std::string(file_out);

    std::string file_name;

    if (type == _util.k_mc){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_mc_run%s.root", _util.run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
            f_nuexsec_tree = new TFile( file_name.c_str(), "UPDATE");
            f_nuexsec_tree->cd();
        }

        // Create the TTree
        if (_util.isfakedata){
            tree             = new TTree("data_tree",                "data_tree");
            nue_tree         = new TTree("data_nue_tree",            "data_nue_tree");
            dedx_tree        = new TTree("data_dedx_tree",           "data_dedx_tree");
            counter_tree     = new TTree("data_counter_tree",        "data_counter_tree");
            nue_counter_tree = new TTree("data_nue_counter_tree",    "data_nue_counter_tree");
        }
        else {
            tree             = new TTree("mc_tree",                "mc_tree");
            nue_tree         = new TTree("mc_nue_tree",            "mc_nue_tree");
            dedx_tree        = new TTree("mc_dedx_tree",           "mc_dedx_tree");
            counter_tree     = new TTree("mc_counter_tree",        "mc_counter_tree");
            nue_counter_tree = new TTree("mc_nue_counter_tree",    "mc_nue_counter_tree");
        }
        
    }
    else if (type == _util.k_data){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_data_run%s.root", _util.run_period);
        else file_name = "files/trees/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec_tree = new TFile(file_name.c_str(), "UPDATE");
            f_nuexsec_tree->cd();
        }

        // Create the TTree
        tree         = new TTree("data_tree",     "data_tree");
        dedx_tree    = new TTree("data_dedx_tree","data_dedx_tree");
        counter_tree = new TTree("data_counter_tree",    "data_counter_tree");

    }
    else if (type == _util.k_ext){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_ext_run%s.root", _util.run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec_tree = new TFile(file_name.c_str(), "UPDATE");
            f_nuexsec_tree->cd();
        }

        // Create the TTree
        tree         = new TTree("ext_tree",     "ext_tree");
        dedx_tree    = new TTree("ext_dedx_tree","ext_dedx_tree");
        counter_tree = new TTree("ext_counter_tree",    "ext_counter_tree");

    }
    else if (type == _util.k_dirt){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_dirt_run%s.root", _util.run_period);
        else file_name = "files/trees/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec_tree = new TFile(file_name.c_str(), "UPDATE");
            f_nuexsec_tree->cd();
        }

        // Create the TTree
        tree         = new TTree("dirt_tree",     "dirt_tree");
        dedx_tree    = new TTree("dirt_dedx_tree","dirt_dedx_tree");
        counter_tree = new TTree("dirt_counter_tree",    "dirt_counter_tree");

    }
    else {
        std::cout << "Unknown input type!! "<<  __PRETTY_FUNCTION__ << std::endl;
        exit(1);
    }

    // Set the type
    _type = type;


    if (std::string(_util.intrinsic_mode) == "intrinsic")
        SetBranches(nue_tree);
    else
        SetBranches(tree);

    dedx_tree->Branch("shr_hits_u_tot", &shr_hits_u_tot, "shr_hits_u_tot/I");
    dedx_tree->Branch("shr_hits_v_tot", &shr_hits_v_tot, "shr_hits_v_tot/I");
    dedx_tree->Branch("shr_hits_y_tot", &shr_hits_y_tot, "shr_hits_y_tot/I");
    
    dedx_tree->Branch("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y, "shr_tkfit_dedx_Y/F");
    dedx_tree->Branch("shr_tkfit_dedx_V", &shr_tkfit_dedx_V, "shr_tkfit_dedx_V/F");
    dedx_tree->Branch("shr_tkfit_dedx_U", &shr_tkfit_dedx_U, "shr_tkfit_dedx_U/F");
    
    dedx_tree->Branch("weight",          &weight,       "weight/D");
    dedx_tree->Branch("classification",  &classification);
    dedx_tree->Branch("shr_distance",    &shr_distance, "shr_distance/F");
    dedx_tree->Branch("shr_theta",       &shr_theta,    "shr_theta/F");
    dedx_tree->Branch("cut", &cut);
    dedx_tree->Branch("npi0", &npi0);
    dedx_tree->Branch("nu_pdg", &nu_pdg);

    if (type == _util.k_mc){
        // Counter Tree
        nue_counter_tree->Branch("count_nue_cc_qe",  &count_nue_cc_qe);
        nue_counter_tree->Branch("count_nue_cc_res", &count_nue_cc_res);
        nue_counter_tree->Branch("count_nue_cc_dis", &count_nue_cc_dis);
        nue_counter_tree->Branch("count_nue_cc_coh", &count_nue_cc_coh);
        nue_counter_tree->Branch("count_nue_cc_mec", &count_nue_cc_mec);
        
        nue_counter_tree->Branch("count_nuebar_cc_qe",  &count_nuebar_cc_qe);
        nue_counter_tree->Branch("count_nuebar_cc_res", &count_nuebar_cc_res);
        nue_counter_tree->Branch("count_nuebar_cc_dis", &count_nuebar_cc_dis);
        nue_counter_tree->Branch("count_nuebar_cc_coh", &count_nuebar_cc_coh);
        nue_counter_tree->Branch("count_nuebar_cc_mec", &count_nuebar_cc_mec);
        
        nue_counter_tree->Branch("count_nue_cc_infv",      &count_nue_cc_infv);
        nue_counter_tree->Branch("count_nuebar_cc_infv",   &count_nuebar_cc_infv);
        counter_tree->Branch("count_nue_cc_incryo",    &count_nue_cc_incryo);
        counter_tree->Branch("count_nuebar_cc_incryo", &count_nuebar_cc_incryo);
        
        counter_tree    ->Branch("count_numu_cc_qe",  &count_numu_cc_qe);
        counter_tree    ->Branch("count_numu_cc_res", &count_numu_cc_res);
        counter_tree    ->Branch("count_numu_cc_dis", &count_numu_cc_dis);
        counter_tree    ->Branch("count_numu_cc_coh", &count_numu_cc_coh);
        counter_tree    ->Branch("count_numu_cc_mec", &count_numu_cc_mec);
        
        counter_tree    ->Branch("count_numubar_cc_qe",  &count_numubar_cc_qe);
        counter_tree    ->Branch("count_numubar_cc_res", &count_numubar_cc_res);
        counter_tree    ->Branch("count_numubar_cc_dis", &count_numubar_cc_dis);
        counter_tree    ->Branch("count_numubar_cc_coh", &count_numubar_cc_coh);
        counter_tree    ->Branch("count_numubar_cc_mec", &count_numubar_cc_mec);
        
        counter_tree    ->Branch("count_numu_cc_infv",      &count_numu_cc_infv);
        counter_tree    ->Branch("count_numubar_cc_infv",   &count_numubar_cc_infv);
        counter_tree    ->Branch("count_numu_cc_incryo",    &count_numu_cc_incryo);
        counter_tree    ->Branch("count_numubar_cc_incryo", &count_numubar_cc_incryo);


        nue_counter_tree->Branch("count_pi0_nue_cc_nopi0",    &count_pi0_nue_cc_nopi0);
        nue_counter_tree->Branch("count_pi0_nue_cc_pi0",      &count_pi0_nue_cc_pi0);
        nue_counter_tree->Branch("count_pi0_nuebar_cc_nopi0", &count_pi0_nuebar_cc_nopi0);
        nue_counter_tree->Branch("count_pi0_nuebar_cc_pi0",   &count_pi0_nuebar_cc_pi0);
        counter_tree    ->Branch("count_pi0_numu_cc_nopi0",   &count_pi0_numu_cc_nopi0);
        counter_tree    ->Branch("count_pi0_numu_cc_pi0",     &count_pi0_numu_cc_pi0);
        counter_tree    ->Branch("count_pi0_nc_nopi0",        &count_pi0_nc_nopi0);
        counter_tree    ->Branch("count_pi0_nc_pi0",          &count_pi0_nc_pi0);
        
        nue_counter_tree->Branch("count_nue_cc",          &count_nue_cc);
        nue_counter_tree->Branch("count_nuebar_cc",       &count_nuebar_cc);
        counter_tree    ->Branch("count_nu_out_fv",       &count_nu_out_fv);
        counter_tree    ->Branch("count_cosmic",          &count_cosmic);
        counter_tree    ->Branch("count_numu_cc",         &count_numu_cc);
        counter_tree    ->Branch("count_numu_cc_pi0",     &count_numu_cc_pi0);
        counter_tree    ->Branch("count_nc",              &count_nc);
        counter_tree    ->Branch("count_nc_pi0",          &count_nc_pi0);
        nue_counter_tree->Branch("count_unmatched",       &count_unmatched);
        nue_counter_tree->Branch("count_unmatched_nue",   &count_unmatched_nue);
        nue_counter_tree->Branch("count_cosmic_nue",      &count_cosmic_nue);
        nue_counter_tree->Branch("count_unmatched_nuebar",&count_unmatched_nuebar);
        nue_counter_tree->Branch("count_cosmic_nuebar",   &count_cosmic_nuebar);
        nue_counter_tree->Branch("count_thr_nue",         &count_thr_nue);
        nue_counter_tree->Branch("count_thr_nuebar",      &count_thr_nuebar);
        counter_tree    ->Branch("count_total_mc",        &count_total_mc);
    }
        
    
    counter_tree    ->Branch("count_data",         &count_data);
    counter_tree    ->Branch("count_ext",          &count_ext);
    counter_tree    ->Branch("count_dirt",         &count_dirt);

    

    std::cout << "Finished Initalising the Tree Helper..."<< std::endl;

}
// -----------------------------------------------------------------------------
void TreeHelper::FillVars(SliceContainer &SC, bool _passed_selection){

    f_nuexsec_tree->cd();

    run                      = SC.run;
    subrun                   = SC.sub;
    event                    = SC.evt;
    gen                      = SC.is_signal;
    passed_selection         = _passed_selection;
    classification           = SC.classification.first;
    weight                   = SC.cv_weight;
    true_energy              = SC.nu_e;
    reco_energy              = SC.reco_e;
    shr_tkfit_dedx_Y         = SC.shr_tkfit_dedx_Y;
    n_showers                = SC.n_showers;
    n_tracks                 = SC.n_tracks;
    shr_theta                = SC.shr_theta;
    shr_phi                  = SC.shr_phi;
    shr_energy_cali          = SC.shr_energy_cali;
    shrmoliereavg            = SC.shrmoliereavg;
    shr_hits_max             = SC.shr_hits_max;
    elec_e                   = SC.elec_e;
    ppfx_cv                  = SC.ppfx_cv;
    weightSplineTimesTune    = SC.weightSplineTimesTune;
    nu_pdg                   = SC.nu_pdg;
    numi_ang                 = SC.nu_angle;
    shr_bkt_pdg              = SC.shr_bkt_pdg;
    npi0                     = SC.npi0;
    pi0_e                    = SC.pi0_e;
    interaction              = SC.interaction;
    effective_angle          = SC.effective_angle;
    cos_effective_angle      = SC.cos_effective_angle;
    true_effective_angle     = SC.true_effective_angle;
    cos_true_effective_angle = std::cos(SC.true_effective_angle * 3.14159 / 180.0);

    shr_bkt_purity       = SC.shr_bkt_purity;
    shr_bkt_completeness = SC.shr_bkt_completeness;
    shr_bkt_E            = SC.shr_bkt_E;
    all_shr_hits         = *SC.all_shr_hits;
    all_shr_energies     = *SC.all_shr_energies;

    
    if (SC.weightsGenie != NULL) weightsGenie           = *SC.weightsGenie; // If these aren't set by default then bad things happen in memory land
    if (SC.weightsReint != NULL) weightsReint           = *SC.weightsReint;
    if (SC.weightsPPFX  != NULL) weightsPPFX            = *SC.weightsPPFX;
    knobRPAup              = SC.knobRPAup;
    knobCCMECup            = SC.knobCCMECup;
    knobAxFFCCQEup         = SC.knobAxFFCCQEup;
    knobVecFFCCQEup        = SC.knobVecFFCCQEup;
    knobDecayAngMECup      = SC.knobDecayAngMECup;
    knobThetaDelta2Npiup   = SC.knobThetaDelta2Npiup;
    knobThetaDelta2NRadup  = SC.knobThetaDelta2NRadup;
    knobRPA_CCQE_Reducedup = SC.knobRPA_CCQE_Reducedup;
    knobNormCCCOHup        = SC.knobNormCCCOHup;
    knobNormNCCOHup        = SC.knobNormNCCOHup;
    knobxsr_scc_Fv3up      = SC.knobxsr_scc_Fv3up;
    knobxsr_scc_Fa3up      = SC.knobxsr_scc_Fa3up;
    knobRPAdn              = SC.knobRPAdn;
    knobCCMECdn            = SC.knobCCMECdn;
    knobAxFFCCQEdn         = SC.knobAxFFCCQEdn;
    knobVecFFCCQEdn        = SC.knobVecFFCCQEdn;
    knobDecayAngMECdn      = SC.knobDecayAngMECdn;
    knobThetaDelta2Npidn   = SC.knobThetaDelta2Npidn;
    knobThetaDelta2NRaddn  = SC.knobThetaDelta2NRaddn;
    knobRPA_CCQE_Reduceddn = SC.knobRPA_CCQE_Reduceddn;
    knobNormCCCOHdn        = SC.knobNormCCCOHdn;
    knobNormNCCOHdn        = SC.knobNormNCCOHdn;
    knobxsr_scc_Fv3dn      = SC.knobxsr_scc_Fv3dn;
    knobxsr_scc_Fa3dn      = SC.knobxsr_scc_Fa3dn;

    // Fill the nue tree
    if (std::string(_util.intrinsic_mode) == "intrinsic")
        nue_tree->Fill();
    else
        tree->Fill();

}
// -----------------------------------------------------------------------------
void TreeHelper::Fill_dedxVars(SliceContainer &SC, std::pair<std::string, int> _classification, std::string _cut, double _weight){

    f_nuexsec_tree->cd();

    classification = _classification.first;
    weight = _weight;
    cut    = _cut;

    npi0 = SC.npi0;

    nu_pdg = SC.nu_pdg;


    shr_tkfit_dedx_Y       = SC.shr_tkfit_dedx_Y;
    shr_tkfit_dedx_V       = SC.shr_tkfit_dedx_V;
    shr_tkfit_dedx_U       = SC.shr_tkfit_dedx_U;

    // Redefine some cases where we have undefined dedx for out study
    if (shr_tkfit_dedx_Y < 0) {
        shr_hits_y_tot = 0;
    }
    else {
        shr_hits_y_tot         = SC.shr_hits_y_tot;
    }

    if (shr_tkfit_dedx_V < 0) {
        shr_hits_v_tot = 0;
    }
    else {
        shr_hits_v_tot         = SC.shr_hits_v_tot;
    }

    if (shr_tkfit_dedx_U < 0) {
        shr_hits_u_tot = 0;
    }
    else {
        shr_hits_u_tot         = SC.shr_hits_u_tot;
    }
    
    
    shr_distance           = SC.shr_distance;
    shr_theta              = SC.shr_theta;

    dedx_tree->Fill();

}
// -----------------------------------------------------------------------------
void TreeHelper::WriteTree(int type){

    f_nuexsec_tree->cd();

    if (std::string(_util.intrinsic_mode) == "intrinsic")
        nue_tree->Write("",TObject::kOverwrite);
    else
        tree->Write("",TObject::kOverwrite);

    if (type == _util.k_mc)
        nue_counter_tree->Write("",TObject::kOverwrite);
    
    if (std::string(_util.intrinsic_mode) != "intrinsic")
        counter_tree->Write("",TObject::kOverwrite);


    dedx_tree->Write("",TObject::kOverwrite);

}
// -----------------------------------------------------------------------------
void TreeHelper::Fill_counters(std::vector<double> counter_v, bool bool_use_mc, bool bool_use_ext, bool bool_use_data, bool bool_use_dirt){

    f_nuexsec_tree->cd();

    if (bool_use_mc){
        count_nue_cc_qe  = counter_v.at(_util.k_count_nue_cc_qe);
        count_nue_cc_res = counter_v.at(_util.k_count_nue_cc_res);
        count_nue_cc_dis = counter_v.at(_util.k_count_nue_cc_dis);
        count_nue_cc_coh = counter_v.at(_util.k_count_nue_cc_coh);
        count_nue_cc_mec = counter_v.at(_util.k_count_nue_cc_mec);
        
        count_nuebar_cc_qe  = counter_v.at(_util.k_count_nuebar_cc_qe);
        count_nuebar_cc_res = counter_v.at(_util.k_count_nuebar_cc_res);
        count_nuebar_cc_dis = counter_v.at(_util.k_count_nuebar_cc_dis);
        count_nuebar_cc_coh = counter_v.at(_util.k_count_nuebar_cc_coh);
        count_nuebar_cc_mec = counter_v.at(_util.k_count_nuebar_cc_mec);
        
        count_nue_cc_infv      = counter_v.at(_util.k_count_nue_cc_infv);
        count_nuebar_cc_infv   = counter_v.at(_util.k_count_nuebar_cc_infv);
        count_nue_cc_incryo    = counter_v.at(_util.k_count_nue_cc_incryo);
        count_nuebar_cc_incryo = counter_v.at(_util.k_count_nuebar_cc_incryo);
        
        count_numu_cc_qe  = counter_v.at(_util.k_count_numu_cc_qe);
        count_numu_cc_res = counter_v.at(_util.k_count_numu_cc_res);
        count_numu_cc_dis = counter_v.at(_util.k_count_numu_cc_dis);
        count_numu_cc_coh = counter_v.at(_util.k_count_numu_cc_coh);
        count_numu_cc_mec = counter_v.at(_util.k_count_numu_cc_mec);
        
        count_numubar_cc_qe  = counter_v.at(_util.k_count_numubar_cc_qe);
        count_numubar_cc_res = counter_v.at(_util.k_count_numubar_cc_res);
        count_numubar_cc_dis = counter_v.at(_util.k_count_numubar_cc_dis);
        count_numubar_cc_coh = counter_v.at(_util.k_count_numubar_cc_coh);
        count_numubar_cc_mec = counter_v.at(_util.k_count_numubar_cc_mec);
        
        count_numu_cc_infv      = counter_v.at(_util.k_count_numu_cc_infv);
        count_numubar_cc_infv   = counter_v.at(_util.k_count_numubar_cc_infv);
        count_numu_cc_incryo    = counter_v.at(_util.k_count_numu_cc_incryo);
        count_numubar_cc_incryo = counter_v.at(_util.k_count_numubar_cc_incryo);

        count_pi0_nue_cc_nopi0    = counter_v.at(_util.k_count_pi0_nue_cc_nopi0);
        count_pi0_nue_cc_pi0      = counter_v.at(_util.k_count_pi0_nue_cc_pi0);
        count_pi0_nuebar_cc_nopi0 = counter_v.at(_util.k_count_pi0_nuebar_cc_nopi0);
        count_pi0_nuebar_cc_pi0   = counter_v.at(_util.k_count_pi0_nuebar_cc_pi0);
        count_pi0_numu_cc_nopi0   = counter_v.at(_util.k_count_pi0_numu_cc_nopi0);
        count_pi0_numu_cc_pi0     = counter_v.at(_util.k_count_pi0_numu_cc_pi0);
        count_pi0_nc_nopi0        = counter_v.at(_util.k_count_pi0_nc_nopi0);
        count_pi0_nc_pi0          = counter_v.at(_util.k_count_pi0_nc_pi0);
        
        count_nue_cc           = counter_v.at(_util.k_count_nue_cc);
        count_nuebar_cc        = counter_v.at(_util.k_count_nuebar_cc);
        count_nu_out_fv        = counter_v.at(_util.k_count_nu_out_fv);
        count_cosmic           = counter_v.at(_util.k_count_cosmic);
        count_numu_cc          = counter_v.at(_util.k_count_numu_cc);
        count_numu_cc_pi0      = counter_v.at(_util.k_count_numu_cc_pi0);
        count_nc               = counter_v.at(_util.k_count_nc);
        count_nc_pi0           = counter_v.at(_util.k_count_nc_pi0);
        count_unmatched        = counter_v.at(_util.k_count_unmatched);
        count_unmatched_nue    = counter_v.at(_util.k_count_unmatched_nue);
        count_cosmic_nue       = counter_v.at(_util.k_count_cosmic_nue);
        count_unmatched_nuebar = counter_v.at(_util.k_count_unmatched_nuebar);
        count_cosmic_nuebar    = counter_v.at(_util.k_count_cosmic_nuebar);
        count_thr_nue          = counter_v.at(_util.k_count_thr_nue);
        count_thr_nuebar       = counter_v.at(_util.k_count_thr_nuebar);
        count_total_mc         = counter_v.at(_util.k_count_total_mc);
    }

    if (bool_use_data)  count_data         = counter_v.at(_util.k_count_data);
    
    if (bool_use_ext)   count_ext          = counter_v.at(_util.k_count_ext);
    
    if (bool_use_dirt)  count_dirt         = counter_v.at(_util.k_count_dirt);

    counter_tree->Fill();

    if (bool_use_mc)
        nue_counter_tree->Fill();

}
// -----------------------------------------------------------------------------
void TreeHelper::SetBranches(TTree * tree){
    
    // Set the tree branches
    tree->Branch("run",                    &run,    "run/I");
    tree->Branch("subrun",                 &subrun, "subrun/I");
    tree->Branch("event",                  &event,  "event/I");
    tree->Branch("gen",                    &gen,    "gen/O");
    tree->Branch("passed_selection",       &passed_selection,    "passed_selection/O");
    tree->Branch("weight",                 &weight, "weight/D");
    tree->Branch("true_energy",            &true_energy, "true_energy/D");
    tree->Branch("reco_energy",            &reco_energy, "reco_energy/D");
    tree->Branch("classification",         &classification);
    tree->Branch("shr_tkfit_dedx_Y",       &shr_tkfit_dedx_Y, "shr_tkfit_dedx_Y/F");
    tree->Branch("n_showers",              &n_showers, "n_showers/F");
    tree->Branch("n_tracks",               &n_tracks,  "n_tracks/F");
    tree->Branch("shr_theta",              &shr_theta, "shr_theta/F");
    tree->Branch("shr_phi",                &shr_phi,   "shr_phi/F");
    tree->Branch("shr_energy_cali",        &shr_energy_cali, "shr_energy_cali/F");
    tree->Branch("shrmoliereavg",          &shrmoliereavg, "shrmoliereavg/F");
    tree->Branch("shr_hits_max",           &shr_hits_max,  "shr_hits_max/F");
    tree->Branch("elec_e",                 &elec_e,  "elec_e/F");
    tree->Branch("ppfx_cv",                &ppfx_cv,  "ppfx_cv/F");
    tree->Branch("weightSplineTimesTune",  &weightSplineTimesTune,  "weightSplineTimesTune/F");
    tree->Branch("nu_pdg",                 &nu_pdg,  "nu_pdg/I");
    tree->Branch("numi_ang",               &numi_ang,  "numi_ang/F");
    tree->Branch("shr_bkt_pdg",            &shr_bkt_pdg,  "shr_bkt_pdg/I");
    tree->Branch("shr_bkt_purity",         &shr_bkt_purity);
    tree->Branch("shr_bkt_completeness",   &shr_bkt_completeness);
    tree->Branch("shr_bkt_E",              &shr_bkt_E);
    tree->Branch("npi0",                   &npi0,  "npi0/I");
    tree->Branch("pi0_e",                  &pi0_e, "pi0_e/D");
    tree->Branch("all_shr_hits",     "std::vector<float>", &all_shr_hits);
    tree->Branch("all_shr_energies", "std::vector<float>", &all_shr_energies);
    tree->Branch("interaction",              &interaction,  "interaction/I");
    tree->Branch("effective_angle",          &effective_angle);
    tree->Branch("cos_effective_angle",      &cos_effective_angle);
    tree->Branch("true_effective_angle",     &true_effective_angle);
    tree->Branch("cos_true_effective_angle", &cos_true_effective_angle);

    tree->Branch("weightsGenie", "std::vector<unsigned short>", &weightsGenie);
    tree->Branch("weightsReint", "std::vector<unsigned short>", &weightsReint);
    tree->Branch("weightsPPFX",  "std::vector<unsigned short>", &weightsPPFX);
    tree->Branch("knobRPAup",             &knobRPAup,             "knobRPAup/D");
    tree->Branch("knobRPAdn",             &knobRPAdn,             "knobRPAdn/D");
    tree->Branch("knobCCMECup",           &knobCCMECup,           "knobCCMECup/D");
    tree->Branch("knobCCMECdn",           &knobCCMECdn,           "knobCCMECdn/D");
    tree->Branch("knobAxFFCCQEup",        &knobAxFFCCQEup,        "knobAxFFCCQEup/D");
    tree->Branch("knobAxFFCCQEdn",        &knobAxFFCCQEdn,        "knobAxFFCCQEdn/D");
    tree->Branch("knobVecFFCCQEup",       &knobVecFFCCQEup,       "knobVecFFCCQEup/D");
    tree->Branch("knobVecFFCCQEdn",       &knobVecFFCCQEdn,       "knobVecFFCCQEdn/D");
    tree->Branch("knobDecayAngMECup",     &knobDecayAngMECup,     "knobDecayAngMECup/D");
    tree->Branch("knobDecayAngMECdn",     &knobDecayAngMECdn,     "knobDecayAngMECdn/D");
    tree->Branch("knobThetaDelta2Npiup",  &knobThetaDelta2Npiup,  "knobThetaDelta2Npiup/D");
    tree->Branch("knobThetaDelta2Npidn",  &knobThetaDelta2Npidn,  "knobThetaDelta2Npidn/D");
    tree->Branch("knobThetaDelta2NRadup", &knobThetaDelta2NRadup, "knobThetaDelta2NRadup/D");
    tree->Branch("knobThetaDelta2NRaddn", &knobThetaDelta2NRaddn, "knobThetaDelta2NRaddn/D");
    tree->Branch("knobRPA_CCQE_Reducedup",&knobRPA_CCQE_Reducedup,"knobRPA_CCQE_Reducedup/D");
    tree->Branch("knobRPA_CCQE_Reduceddn",&knobRPA_CCQE_Reduceddn,"knobRPA_CCQE_Reduceddn/D");
    tree->Branch("knobNormCCCOHup",       &knobNormCCCOHup,       "knobNormCCCOHup/D");
    tree->Branch("knobNormCCCOHdn",       &knobNormCCCOHdn,       "knobNormCCCOHdn/D");
    tree->Branch("knobNormNCCOHup",       &knobNormNCCOHup,       "knobNormNCCOHup/D");
    tree->Branch("knobNormNCCOHdn",       &knobNormNCCOHdn,       "knobNormNCCOHdn/D");
    tree->Branch("knobxsr_scc_Fv3up",     &knobxsr_scc_Fv3up,     "knobxsr_scc_Fv3up/D");
    tree->Branch("knobxsr_scc_Fv3dn",     &knobxsr_scc_Fv3dn,     "knobxsr_scc_Fv3dn/D");
    tree->Branch("knobxsr_scc_Fa3up",     &knobxsr_scc_Fa3up,     "knobxsr_scc_Fa3up/D");
    tree->Branch("knobxsr_scc_Fa3dn",     &knobxsr_scc_Fa3dn,     "knobxsr_scc_Fa3dn/D");

}
// -----------------------------------------------------------------------------
