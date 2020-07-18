#include "../include/TreeHelper.h"

// -----------------------------------------------------------------------------
void TreeHelper::Initialise(int type, const char* run_period, const char * file_out ){

    std::cout << "Initalising Tree Helper..." << std::endl;

    std::string file_out_str = file_out;

    std::string file_name;

    if (type == _util.k_mc){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_mc_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
            f_nuexsec = new TFile( file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("mc_tree",         "mc_tree");
        dedx_tree = new TTree("mc_dedx_tree",    "mc_dedx_tree");
        counter_tree = new TTree("mc_counter_tree",    "mc_counter_tree");
    }
    else if (type == _util.k_data){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_data_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("data_tree",     "data_tree");
        dedx_tree = new TTree("data_dedx_tree","data_dedx_tree");
        counter_tree = new TTree("data_counter_tree",    "data_counter_tree");

    }
    else if (type == _util.k_ext){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_ext_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("ext_tree",     "ext_tree");
        dedx_tree = new TTree("ext_dedx_tree","ext_dedx_tree");
        counter_tree = new TTree("ext_counter_tree",    "ext_counter_tree");

    }
    else if (type == _util.k_dirt){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_dirt_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("dirt_tree",     "dirt_tree");
        dedx_tree = new TTree("dirt_dedx_tree","dirt_dedx_tree");
        counter_tree = new TTree("dirt_counter_tree",    "dirt_counter_tree");

    }
    else {
        std::cout << "Unknown input type!! "<<  __PRETTY_FUNCTION__ << std::endl;
        exit(1);
    }

    // Set the type
    _type = type;
    
    // Set the tree branches
    tree->Branch("run",    &run,    "run/I");
    tree->Branch("subrun", &subrun, "subrun/I");
    tree->Branch("event",  &event,  "event/I");
    tree->Branch("gen",    &gen,    "gen/O");
    tree->Branch("weight", &weight, "weight/D");
    tree->Branch("true_energy", &true_energy, "true_energy/D");
    tree->Branch("reco_energy", &reco_energy, "reco_energy/D");
    tree->Branch("classifcation",   &classifcation);
    tree->Branch("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y, "shr_tkfit_dedx_Y/F");
    tree->Branch("n_showers", &n_showers, "n_showers/F");
    tree->Branch("n_tracks",  &n_tracks,  "n_tracks/F");
    tree->Branch("shr_theta", &shr_theta, "shr_theta/F");
    tree->Branch("shr_phi",   &shr_phi,   "shr_phi/F");
    tree->Branch("shr_energy_tot_cali", &shr_energy_tot_cali, "shr_energy_tot_cali/F");
    tree->Branch("shrmoliereavg", &shrmoliereavg, "shrmoliereavg/F");
    tree->Branch("shr_hits_max",  &shr_hits_max,  "shr_hits_max/F");
    tree->Branch("elec_e",  &elec_e,  "elec_e/F");

    tree->Branch("weightsGenie", "std::vector<unsigned short>", &weightsGenie);
    tree->Branch("weightsReint", "std::vector<unsigned short>", &weightsReint);
    tree->Branch("weightsPPFX", "std::vector<unsigned short>", &weightsPPFX);
    tree->Branch("knobRPAup",&knobRPAup,"knobRPAup/D");
    tree->Branch("knobRPAdn",&knobRPAdn,"knobRPAdn/D");
    tree->Branch("knobCCMECup",&knobCCMECup,"knobCCMECup/D");
    tree->Branch("knobCCMECdn",&knobCCMECdn,"knobCCMECdn/D");
    tree->Branch("knobAxFFCCQEup",&knobAxFFCCQEup,"knobAxFFCCQEup/D");
    tree->Branch("knobAxFFCCQEdn",&knobAxFFCCQEdn,"knobAxFFCCQEdn/D");
    tree->Branch("knobVecFFCCQEup",&knobVecFFCCQEup,"knobVecFFCCQEup/D");
    tree->Branch("knobVecFFCCQEdn",&knobVecFFCCQEdn,"knobVecFFCCQEdn/D");
    tree->Branch("knobDecayAngMECup",&knobDecayAngMECup,"knobDecayAngMECup/D");
    tree->Branch("knobDecayAngMECdn",&knobDecayAngMECdn,"knobDecayAngMECdn/D");
    tree->Branch("knobThetaDelta2Npiup",&knobThetaDelta2Npiup,"knobThetaDelta2Npiup/D");
    tree->Branch("knobThetaDelta2Npidn",&knobThetaDelta2Npidn,"knobThetaDelta2Npidn/D");
    tree->Branch("knobThetaDelta2NRadup",&knobThetaDelta2NRadup,"knobThetaDelta2NRadup/D");
    tree->Branch("knobThetaDelta2NRaddn",&knobThetaDelta2NRaddn,"knobThetaDelta2NRaddn/D");
    tree->Branch("knobRPA_CCQE_Reducedup",&knobRPA_CCQE_Reducedup,"knobRPA_CCQE_Reducedup/D");
    tree->Branch("knobRPA_CCQE_Reduceddn",&knobRPA_CCQE_Reduceddn,"knobRPA_CCQE_Reduceddn/D");
    tree->Branch("knobNormCCCOHup",&knobNormCCCOHup,"knobNormCCCOHup/D");
    tree->Branch("knobNormCCCOHdn",&knobNormCCCOHdn,"knobNormCCCOHdn/D");
    tree->Branch("knobNormNCCOHup",&knobNormNCCOHup,"knobNormNCCOHup/D");
    tree->Branch("knobNormNCCOHdn",&knobNormNCCOHdn,"knobNormNCCOHdn/D");

    dedx_tree->Branch("shr_dedx_Y_cali", &shr_dedx_Y_cali, "shr_dedx_Y_cali/F");
    dedx_tree->Branch("shr_dedx_V_cali", &shr_dedx_V_cali, "shr_dedx_V_cali/F");
    dedx_tree->Branch("shr_dedx_U_cali", &shr_dedx_U_cali, "shr_dedx_U_cali/F");
    
    dedx_tree->Branch("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y, "shr_tkfit_dedx_Y/F");
    dedx_tree->Branch("shr_tkfit_dedx_V", &shr_tkfit_dedx_V, "shr_tkfit_dedx_V/F");
    dedx_tree->Branch("shr_tkfit_dedx_U", &shr_tkfit_dedx_U, "shr_tkfit_dedx_U/F");
    
    dedx_tree->Branch("shr_tkfit_dedx_Y_alt", &shr_tkfit_dedx_Y_alt, "shr_tkfit_dedx_Y_alt/F");
    dedx_tree->Branch("shr_tkfit_dedx_V_alt", &shr_tkfit_dedx_V_alt, "shr_tkfit_dedx_V_alt/F");
    dedx_tree->Branch("shr_tkfit_dedx_U_alt", &shr_tkfit_dedx_U_alt, "shr_tkfit_dedx_U_alt/F");
    
    dedx_tree->Branch("shr_tkfit_2cm_dedx_Y", &shr_tkfit_2cm_dedx_Y, "shr_tkfit_2cm_dedx_Y/F");
    dedx_tree->Branch("shr_tkfit_2cm_dedx_V", &shr_tkfit_2cm_dedx_V, "shr_tkfit_2cm_dedx_V/F");
    dedx_tree->Branch("shr_tkfit_2cm_dedx_U", &shr_tkfit_2cm_dedx_U, "shr_tkfit_2cm_dedx_U/F");
    
    dedx_tree->Branch("shr_tkfit_gap05_dedx_Y", &shr_tkfit_gap05_dedx_Y, "shr_tkfit_gap05_dedx_Y/F");
    dedx_tree->Branch("shr_tkfit_gap05_dedx_V", &shr_tkfit_gap05_dedx_V, "shr_tkfit_gap05_dedx_V/F");
    dedx_tree->Branch("shr_tkfit_gap05_dedx_U", &shr_tkfit_gap05_dedx_U, "shr_tkfit_gap05_dedx_U/F");
    
    dedx_tree->Branch("shr_tkfit_gap10_dedx_Y", &shr_tkfit_gap10_dedx_Y, "shr_tkfit_gap10_dedx_Y/F");
    dedx_tree->Branch("shr_tkfit_gap10_dedx_V", &shr_tkfit_gap10_dedx_V, "shr_tkfit_gap10_dedx_V/F");
    dedx_tree->Branch("shr_tkfit_gap10_dedx_U", &shr_tkfit_gap10_dedx_U, "shr_tkfit_gap10_dedx_U/F");

    dedx_tree->Branch("weight",          &weight,       "weight/D");
    dedx_tree->Branch("classifcation",   &classifcation);
    dedx_tree->Branch("shr_distance",    &shr_distance, "shr_distance/F");
    dedx_tree->Branch("shr_theta",       &shr_theta,    "shr_theta/F");
    dedx_tree->Branch("cut", &cut);


    // Counter Tree
    counter_tree->Branch("count_nue_cc_qe",  &count_nue_cc_qe);
    counter_tree->Branch("count_nue_cc_res", &count_nue_cc_res);
    counter_tree->Branch("count_nue_cc_dis", &count_nue_cc_dis);
    counter_tree->Branch("count_nue_cc_coh", &count_nue_cc_coh);
    counter_tree->Branch("count_nue_cc_mec", &count_nue_cc_mec);
    
    counter_tree->Branch("count_nuebar_cc_qe",  &count_nuebar_cc_qe);
    counter_tree->Branch("count_nuebar_cc_res", &count_nuebar_cc_res);
    counter_tree->Branch("count_nuebar_cc_dis", &count_nuebar_cc_dis);
    counter_tree->Branch("count_nuebar_cc_coh", &count_nuebar_cc_coh);
    counter_tree->Branch("count_nuebar_cc_mec", &count_nuebar_cc_mec);
    
    counter_tree->Branch("count_nue_cc_infv",      &count_nue_cc_infv);
    counter_tree->Branch("count_nuebar_cc_infv",   &count_nuebar_cc_infv);
    counter_tree->Branch("count_nue_cc_incryo",    &count_nue_cc_incryo);
    counter_tree->Branch("count_nuebar_cc_incryo", &count_nuebar_cc_incryo);
    
    counter_tree->Branch("count_numu_cc_qe",  &count_numu_cc_qe);
    counter_tree->Branch("count_numu_cc_res", &count_numu_cc_res);
    counter_tree->Branch("count_numu_cc_dis", &count_numu_cc_dis);
    counter_tree->Branch("count_numu_cc_coh", &count_numu_cc_coh);
    counter_tree->Branch("count_numu_cc_mec", &count_numu_cc_mec);
    
    counter_tree->Branch("count_numubar_cc_qe",  &count_numubar_cc_qe);
    counter_tree->Branch("count_numubar_cc_res", &count_numubar_cc_res);
    counter_tree->Branch("count_numubar_cc_dis", &count_numubar_cc_dis);
    counter_tree->Branch("count_numubar_cc_coh", &count_numubar_cc_coh);
    counter_tree->Branch("count_numubar_cc_mec", &count_numubar_cc_mec);
    
    counter_tree->Branch("count_numu_cc_infv",      &count_numu_cc_infv);
    counter_tree->Branch("count_numubar_cc_infv",   &count_numubar_cc_infv);
    counter_tree->Branch("count_numu_cc_incryo",    &count_numu_cc_incryo);
    counter_tree->Branch("count_numubar_cc_incryo", &count_numubar_cc_incryo);
    
    counter_tree->Branch("count_nue_cc",       &count_nue_cc);
    counter_tree->Branch("count_nuebar_cc",    &count_nuebar_cc);
    counter_tree->Branch("count_nu_out_fv",    &count_nu_out_fv);
    counter_tree->Branch("count_cosmic",       &count_cosmic);
    counter_tree->Branch("count_numu_cc",      &count_numu_cc);
    counter_tree->Branch("count_numu_cc_pi0",  &count_numu_cc_pi0);
    counter_tree->Branch("count_nc",           &count_nc);
    counter_tree->Branch("count_nc_pi0",       &count_nc_pi0);
    counter_tree->Branch("count_unmatched",    &count_unmatched);
    counter_tree->Branch("count_total_mc",     &count_total_mc);
    counter_tree->Branch("count_data",         &count_data);
    counter_tree->Branch("count_ext",          &count_ext);
    counter_tree->Branch("count_dirt",         &count_dirt);

    std::cout << "Finished Initalising the Tree Helper..."<< std::endl;

}
// -----------------------------------------------------------------------------
void TreeHelper::FillVars(SliceContainer &SC, std::pair<std::string, int> _classification, bool _gen, double _weight, double _reco_energy){

    f_nuexsec->cd();

    run    = SC.run;
    subrun = SC.sub;
    event  = SC.evt;
    gen    = _gen;
    classifcation = _classification.first;
    weight = _weight;
    true_energy = SC.nu_e;
    reco_energy = _reco_energy;
    shr_tkfit_dedx_Y = SC.shr_tkfit_dedx_Y;
    n_showers = SC.n_showers;
    n_tracks  = SC.n_tracks;
    shr_theta = SC.shr_theta;
    shr_phi   = SC.shr_phi;
    shr_energy_tot_cali = SC.shr_energy_tot_cali;
    shrmoliereavg = SC.shrmoliereavg;
    shr_hits_max  = SC.shr_hits_max;
    elec_e   = SC.elec_e;

    
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


    tree->Fill();

}
// -----------------------------------------------------------------------------
void TreeHelper::Fill_dedxVars(SliceContainer &SC, std::pair<std::string, int> _classification, std::string _cut, double _weight){

    f_nuexsec->cd();

    classifcation = _classification.first;
    weight = _weight;
    cut    = _cut;


    shr_dedx_Y_cali        = SC.shr_dedx_Y_cali;
    shr_dedx_V_cali        = SC.shr_dedx_V_cali;
    shr_dedx_U_cali        = SC.shr_dedx_U_cali;
    shr_tkfit_dedx_Y       = SC.shr_tkfit_dedx_Y;
    shr_tkfit_dedx_V       = SC.shr_tkfit_dedx_V;
    shr_tkfit_dedx_U       = SC.shr_tkfit_dedx_U;
    shr_tkfit_dedx_Y_alt   = SC.shr_tkfit_dedx_Y_alt;
    shr_tkfit_dedx_V_alt   = SC.shr_tkfit_dedx_V_alt;
    shr_tkfit_dedx_U_alt   = SC.shr_tkfit_dedx_U_alt;
    shr_tkfit_2cm_dedx_Y   = SC.shr_tkfit_2cm_dedx_Y;
    shr_tkfit_2cm_dedx_V   = SC.shr_tkfit_2cm_dedx_V;
    shr_tkfit_2cm_dedx_U   = SC.shr_tkfit_2cm_dedx_U;
    shr_tkfit_gap05_dedx_Y = SC.shr_tkfit_gap05_dedx_Y;
    shr_tkfit_gap05_dedx_V = SC.shr_tkfit_gap05_dedx_V;
    shr_tkfit_gap05_dedx_U = SC.shr_tkfit_gap05_dedx_U;
    shr_tkfit_gap10_dedx_Y = SC.shr_tkfit_gap10_dedx_Y;
    shr_tkfit_gap10_dedx_V = SC.shr_tkfit_gap10_dedx_V;
    shr_tkfit_gap10_dedx_U = SC.shr_tkfit_gap10_dedx_U;
    shr_distance           = SC.shr_distance;
    shr_theta              = SC.shr_theta;

    dedx_tree->Fill();

}
// -----------------------------------------------------------------------------
void TreeHelper::WriteTree(){

    f_nuexsec->cd();

    tree->Write("",TObject::kOverwrite);

    dedx_tree->Write("",TObject::kOverwrite);

    counter_tree->Write("",TObject::kOverwrite);

}
// -----------------------------------------------------------------------------
void TreeHelper::Fill_counters(std::vector<double> counter_v, bool bool_use_mc, bool bool_use_ext, bool bool_use_data, bool bool_use_dirt){

    f_nuexsec->cd();

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
        
        count_nue_cc       = counter_v.at(_util.k_count_nue_cc);
        count_nuebar_cc    = counter_v.at(_util.k_count_nuebar_cc);
        count_nu_out_fv    = counter_v.at(_util.k_count_nu_out_fv);
        count_cosmic       = counter_v.at(_util.k_count_cosmic);
        count_numu_cc      = counter_v.at(_util.k_count_numu_cc);
        count_numu_cc_pi0  = counter_v.at(_util.k_count_numu_cc_pi0);
        count_nc           = counter_v.at(_util.k_count_nc);
        count_nc_pi0       = counter_v.at(_util.k_count_nc_pi0);
        count_unmatched    = counter_v.at(_util.k_count_unmatched);
        count_total_mc     = counter_v.at(_util.k_count_total_mc);
    }

    if (bool_use_data)  count_data         = counter_v.at(_util.k_count_data);
    
    if (bool_use_ext)   count_ext          = counter_v.at(_util.k_count_ext);
    
    if (bool_use_dirt)  count_dirt         = counter_v.at(_util.k_count_dirt);

    counter_tree->Fill();

}