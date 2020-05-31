#include "../include/TreeHelper.h"

// -----------------------------------------------------------------------------
void TreeHelper::Initialise(int type, const char* run_period, const char * file_out, int weight_cfg ){

    std::cout << "Initalising Tree Helper..." << std::endl;

    std::string file_out_str = file_out;

    std::string file_name;

    if (type == _util.k_mc){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/mcc8_nuexsec_selected_tree_mc_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
            f_nuexsec = new TFile( file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("mc_tree",         "mc_tree");
        counter_tree = new TTree("mc_counter_tree",    "mc_counter_tree");
        tree->SetDirectory(0);
        counter_tree->SetDirectory(0);
    }
    else if (type == _util.k_data){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/mcc8_nuexsec_selected_tree_data_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("data_tree",     "data_tree");
        counter_tree = new TTree("data_counter_tree",    "data_counter_tree");

    }
    else if (type == _util.k_ext){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/mcc8_nuexsec_selected_tree_ext_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("ext_tree",     "ext_tree");
        counter_tree = new TTree("ext_counter_tree",    "ext_counter_tree");

    }
    else if (type == _util.k_dirt){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/trees/mcc8_nuexsec_selected_tree_dirt_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

        // Create the TTree
        tree      = new TTree("dirt_tree",     "dirt_tree");
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
    tree->Branch("classifcation",   &classifcation);



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
    counter_tree->Branch("count_nue_cc_mixed", &count_nue_cc_mixed);
    counter_tree->Branch("count_nue_cc_out_fv",&count_nue_cc_out_fv);
    counter_tree->Branch("count_cosmic",       &count_cosmic);
    counter_tree->Branch("count_numu_cc",      &count_numu_cc);
    counter_tree->Branch("count_nc_mixed",     &count_nc_mixed);
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
void TreeHelper::FillVars(SliceContainer &SC, std::pair<std::string, int> _classification, bool _gen){

    f_nuexsec->cd();

    run    = SC.run;
    subrun = SC.subrun;
    event  = SC.event;
    gen    = _gen;
    classifcation = _classification.first;
   
    tree->Fill();

}
// -----------------------------------------------------------------------------
void TreeHelper::WriteTree(int type){

    f_nuexsec->cd();

    tree->Write("",TObject::kOverwrite);

    counter_tree->Write("",TObject::kOverwrite);
    

}
// -----------------------------------------------------------------------------
void TreeHelper::Fill_counters(std::vector<double> counter_v, std::string cut_name, bool bool_use_mc, bool bool_use_ext, bool bool_use_data, bool bool_use_dirt){

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
        count_nue_cc_mixed = counter_v.at(_util.k_count_nue_cc_mixed);
        count_nue_cc_out_fv = counter_v.at(_util.k_count_nue_cc_out_fv);
        count_cosmic       = counter_v.at(_util.k_count_cosmic);
        count_numu_cc      = counter_v.at(_util.k_count_numu_cc);
        count_nc_mixed     = counter_v.at(_util.k_count_nc_mixed);
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