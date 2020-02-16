#include "../include/utility.h"


// -----------------------------------------------------------------------------
std::vector<double> utility::configure_cuts(
        double _x1,
        double _x2,
        double _y1,
        double _y2,
        double _z1,
        double _z2,
        double flash_pe_threshold,
        double flash_time_start,
        double flash_time_end,
        double tolerance,
        double shwr_nue_tolerance,
        double trk_nue_tolerance,
        double shwr_hit_threshold,
        double shwr_hit_threshold_collection,
        double tolerance_open_angle_min,
        double tolerance_open_angle_max,
        double tolerance_dedx_min,
        double tolerance_dedx_max,
        double dist_tolerance,
        double pfp_hits_length_tolerance,
        double ratio_tolerance,
        bool   do_variations
        ) {
    std::vector<double> config;
    config.resize(22,0);

    config[0]  = _x1;
    config[1]  = _x2;
    config[2]  = _y1;
    config[3]  = _y2;
    config[4]  = _z1;
    config[5]  = _z2;
    config[6]  = flash_pe_threshold;
    config[7]  = flash_time_start;
    config[8]  = flash_time_end;
    config[9]  = tolerance;
    config[10] = shwr_nue_tolerance;
    config[11] = trk_nue_tolerance;
    config[12] = shwr_hit_threshold;
    config[13] = shwr_hit_threshold_collection;
    config[14] = tolerance_open_angle_min;
    config[15] = tolerance_open_angle_max;
    config[16] = tolerance_dedx_min;
    config[17] = tolerance_dedx_max;
    config[18] = dist_tolerance;
    config[19] = pfp_hits_length_tolerance;
    config[20] = ratio_tolerance;
    config[21] = do_variations;

    return config;

} // End config function
// -----------------------------------------------------------------------------
bool utility::GetFile(TFile* &f, TString string){
    f = TFile::Open(string);
    
    if (f == NULL) {
        std::cout << "failed to get:\t" << string << "\tThis mode wont be used in the selection" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
void utility::GetTree(TFile* f, TTree* &T, TString string){
    T = (TTree*) f->Get(string);
    if (T == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis tree might not exist in the file, exiting...\n" << std::endl;
        exit(1);
    }
    else {
        return;
    }
}
// -----------------------------------------------------------------------------
bool utility::GetHist(TFile* f, TH1D* &h, TString string){
    h = (TH1D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool utility::GetDirectory(TFile* f, TDirectory* &d, TString string){
    d = (TDirectory*)f->Get(string);
    if (d == NULL) {
        // std::cout << "\nfailed to get:\t" << string << "\tThis directory might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
utility::utility(){
    // For creating histogram names
    std::vector<std::string> type_prefix = {"MC", "Data", "EXT", "Dirt"};

    // Cut directory names
    std::vector<std::string> cut_dirs = {
            "In_FV"                            // Fiducial volume
            };


    // Names of the plot types
    std::vector<std::string> plot_types = {
                "Truth",
                "Reco",
                "Optical",
                "Stack"
                };

    // Names of the classifications
    std::vector<std::string> classification_dirs = {
                "nue_cc",
                "nue_cc_mixed",
                "nu_out_fv",
                "cosmic",
                "numu_cc",
                "numu_cc_pi0",
                "nc",
                "nc_pi0",
                "unmatched",
                "ext",
                "data",
                "dirt"
                };

}
// -----------------------------------------------------------------------------
void utility::Tabulate(std::string interaction, std::string classification, int type, std::vector<int> &counter_v) {
    
    if (type == k_mc){
        
        // Define the counters

        // nue interaction types
        int nue_cc_qe     = 0;
        int nue_cc_res    = 0;
        int nue_cc_dis    = 0;
        int nue_cc_coh    = 0;
        int nue_cc_mec    = 0;

        // nue bar interaction types
        int nue_bar_cc_qe     = 0;
        int nue_bar_cc_res    = 0;
        int nue_bar_cc_dis    = 0;
        int nue_bar_cc_coh    = 0;
        int nue_bar_cc_mec    = 0;

        // total nue + nuebar interaction types
        int total_nue_cc_qe     = 0;
        int total_nue_cc_res    = 0;
        int total_nue_cc_dis    = 0;
        int total_nue_cc_coh    = 0;
        int total_nue_cc_mec    = 0;

        // total numu + numubar interaction types
        int numu_cc_qe    = 0;
        int numu_cc_res   = 0;
        int numu_cc_dis   = 0;
        int numu_cc_coh   = 0;
        int numu_cc_mec   = 0;
        
        // Selection Categories
        int nue_cc        = 0;
        int nue_cc_mixed  = 0;
        int cosmic        = 0;
        int numu_cc       = 0;
        int numu_cc_pi0   = 0;
        int nu_out_fv     = 0;
        int nc            = 0;
        int nc_pi0        = 0;
        int unmatched     = 0;
        int total         = 0;

        // Total nue or nuebar interaction types
        int only_nue_cc     = 0;
        int only_nue_bar_cc = 0;

        if (interaction == "nue_cc_qe")  nue_cc_qe++; 
        if (interaction == "nue_cc_res") nue_cc_res++;
        if (interaction == "nue_cc_coh") nue_cc_coh++;
        if (interaction == "nue_cc_dis") nue_cc_dis++;
        if (interaction == "nue_cc_mec") nue_cc_mec++;

        if (interaction == "nue_bar_cc_qe")  nue_bar_cc_qe++; 
        if (interaction == "nue_bar_cc_res") nue_bar_cc_res++;
        if (interaction == "nue_bar_cc_coh") nue_bar_cc_coh++;
        if (interaction == "nue_bar_cc_dis") nue_bar_cc_dis++;
        if (interaction == "nue_bar_cc_mec") nue_bar_cc_mec++;

        if (classification == "numu_cc_qe"  || classification == "numu_bar_cc_qe")  numu_cc_qe++;
        if (classification == "numu_cc_res" || classification == "numu_bar_cc_res") numu_cc_res++;
        if (classification == "numu_cc_coh" || classification == "numu_bar_cc_coh") numu_cc_coh++;
        if (classification == "numu_cc_dis" || classification == "numu_bar_cc_dis") numu_cc_dis++;
        if (classification == "numu_cc_mec" || classification == "numu_bar_cc_mec") numu_cc_mec++;
        if (classification == "numu_cc_pi0" || classification == "numu_bar_cc_pi0") numu_cc_pi0++;
        
        if (classification == "nu_out_fv")    nu_out_fv++;
        if (classification == "nue_cc_mixed") nue_cc_mixed++;
        if (classification == "nc")           nc++;
        if (classification == "nc_pi0")       nc_pi0++;
        if (classification == "cosmic")       cosmic++;
        if (classification == "unmatched")    unmatched++;

        only_nue_cc     = nue_cc_qe     + nue_cc_res     + nue_cc_dis     + nue_cc_coh     + nue_cc_mec;
        only_nue_bar_cc = nue_bar_cc_qe + nue_bar_cc_res + nue_bar_cc_dis + nue_bar_cc_coh + nue_bar_cc_mec;

        total_nue_cc_qe  = nue_cc_qe  + nue_bar_cc_qe;
        total_nue_cc_res = nue_cc_res + nue_bar_cc_res;
        total_nue_cc_dis = nue_cc_dis + nue_bar_cc_dis;
        total_nue_cc_coh = nue_cc_coh + nue_bar_cc_coh;
        total_nue_cc_mec = nue_cc_mec + nue_bar_cc_mec;

        // Total nue + nuebar
        nue_cc = only_nue_cc + only_nue_bar_cc;

        // Total numu's
        numu_cc = numu_cc_qe + numu_cc_res + numu_cc_dis + numu_cc_coh + numu_cc_mec;
        
        // Add up all the MC events selected including the backgrounds
        total = nue_cc + nue_cc_mixed + nu_out_fv + cosmic + numu_cc + numu_cc_pi0 + nc + nc_pi0 + unmatched;
        
        // Now add counters to the vector
        counter_v.at(k_count_total_nue_cc_qe)  += total_nue_cc_qe;
        counter_v.at(k_count_total_nue_cc_res) += total_nue_cc_res;
        counter_v.at(k_count_total_nue_cc_dis) += total_nue_cc_dis;
        counter_v.at(k_count_total_nue_cc_coh) += total_nue_cc_coh;
        counter_v.at(k_count_total_nue_cc_mec) += total_nue_cc_mec;

        counter_v.at(k_count_only_nue_cc)      += only_nue_cc;
        counter_v.at(k_count_only_nue_bar_cc)  += only_nue_bar_cc;

        counter_v.at(k_count_numu_cc_qe)       += numu_cc_qe;
        counter_v.at(k_count_numu_cc_res)      += numu_cc_res;
        counter_v.at(k_count_numu_cc_dis)      += numu_cc_dis;
        counter_v.at(k_count_numu_cc_coh)      += numu_cc_coh;
        counter_v.at(k_count_numu_cc_mec)      += numu_cc_mec;
        
        counter_v.at(k_count_nue_cc)           += nue_cc;      
        counter_v.at(k_count_nue_cc_mixed)     += nue_cc_mixed;
        counter_v.at(k_count_nu_out_fv)        += nu_out_fv;
        counter_v.at(k_count_cosmic)           += cosmic;
        counter_v.at(k_count_numu_cc)          += numu_cc;
        counter_v.at(k_count_numu_cc_pi0)      += numu_cc_pi0;
        counter_v.at(k_count_nc)               += nc;
        counter_v.at(k_count_nc_pi0)           += nc_pi0;
        counter_v.at(k_count_unmatched)        += unmatched;
        
        counter_v.at(k_count_total)            += total;
    
    }
    else if (type == k_data) {
        counter_v.at(k_count_data) += 1;
    }
    else if (type == k_ext){
        counter_v.at(k_count_ext)  += 1;
    }
    else if (type == k_dirt){
        counter_v.at(k_count_dirt) += 1;
    }
    else {

        std::cout << "unkown type specified!!!  " << __PRETTY_FUNCTION__ << std::endl;
    }
    
}
// -----------------------------------------------------------------------------