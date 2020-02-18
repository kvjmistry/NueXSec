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
utility::utility(){}
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

        int tot_nue_numu_nc = 0;
        
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

        if (interaction == "numu_cc_qe"  || interaction == "numu_bar_cc_qe")  numu_cc_qe++;
        if (interaction == "numu_cc_res" || interaction == "numu_bar_cc_res") numu_cc_res++;
        if (interaction == "numu_cc_coh" || interaction == "numu_bar_cc_coh") numu_cc_coh++;
        if (interaction == "numu_cc_dis" || interaction == "numu_bar_cc_dis") numu_cc_dis++;
        if (interaction == "numu_cc_mec" || interaction == "numu_bar_cc_mec") numu_cc_mec++;
        
        // Check if string contains "NC"
        if (interaction.find("nc") != std::string::npos) {
            tot_nue_numu_nc++;
        }

        if (classification == "nue_cc")       nue_cc++;
        if (classification == "nue_cc_mixed") nue_cc_mixed++;
        if (classification == "nu_out_fv")    nu_out_fv++;
        if (classification == "nc")           nc++;
        if (classification == "nc_pi0")       nc_pi0++;
        if (classification == "numu_cc")      numu_cc++;
        if (classification == "numu_cc_pi0")  numu_cc_pi0++;
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
        // nue_cc = only_nue_cc + only_nue_bar_cc;

        // Total numu's
        // numu_cc = numu_cc_qe + numu_cc_res + numu_cc_dis + numu_cc_coh + numu_cc_mec;
        
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

        counter_v.at(k_count_tot_nue_numu_nc)  += tot_nue_numu_nc;
        
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
void utility::PrintInfo(std::vector<int> counter_v, double intime_scale_factor, double data_scale_factor, double dirt_scale_factor, std::string cut_name, int tot_true_infv_nues) {

    int counter_nue_cc_qe      = counter_v.at(k_count_total_nue_cc_qe);
    int counter_nue_cc_res     = counter_v.at(k_count_total_nue_cc_res);
    int counter_nue_cc_dis     = counter_v.at(k_count_total_nue_cc_dis);
    int counter_nue_cc_coh     = counter_v.at(k_count_total_nue_cc_coh);
    int counter_nue_cc_mec     = counter_v.at(k_count_total_nue_cc_mec);
    
    int counter_numu_cc_qe     = counter_v.at(k_count_numu_cc_qe);
    int counter_numu_cc_res    = counter_v.at(k_count_numu_cc_res);
    int counter_numu_cc_dis    = counter_v.at(k_count_numu_cc_dis);
    int counter_numu_cc_coh    = counter_v.at(k_count_numu_cc_coh);
    int counter_numu_cc_mec    = counter_v.at(k_count_numu_cc_mec);

    int counter_tot_nue_numu_nc = counter_v.at(k_count_tot_nue_numu_nc);

    int counter_nue_cc         = counter_v.at(k_count_nue_cc);
    int counter_nue_cc_mixed   = counter_v.at(k_count_nue_cc_mixed);
    int counter_nu_out_fv      = counter_v.at(k_count_nu_out_fv);
    int counter_cosmic         = counter_v.at(k_count_cosmic);
    int counter_numu_cc        = counter_v.at(k_count_numu_cc);
    int counter_numu_cc_pi0    = counter_v.at(k_count_numu_cc_pi0);
    int counter_nc             = counter_v.at(k_count_nc);
    int counter_nc_pi0         = counter_v.at(k_count_nc_pi0);
    int counter_unmatched      = counter_v.at(k_count_unmatched);
    int counter                = counter_v.at(k_count_total);
    
    int counter_data           = counter_v.at(k_count_data);
    int counter_ext            = counter_v.at(k_count_ext);
    int counter_dirt           = counter_v.at(k_count_dirt);

    counter = counter + (counter_ext * (intime_scale_factor / data_scale_factor)) + (counter_dirt * (dirt_scale_factor / data_scale_factor));

    std::cout << "\n------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "\n\033[0;33m <" << cut_name << "> \033[0m" << std::endl;
    std::cout << " Total Candidate Nue     : " << counter                << "\t \t " << double(counter                * data_scale_factor  ) << std::endl;
    std::cout << " Number of Nue CC        : " << counter_nue_cc         << "\t \t " << double(counter_nue_cc         * data_scale_factor  ) << std::endl;
    std::cout << " Number of Nue CC Mixed  : " << counter_nue_cc_mixed   << "\t \t " << double(counter_nue_cc_mixed   * data_scale_factor  ) << std::endl;
    std::cout << " Number of Nu out FV     : " << counter_nu_out_fv      << "\t \t " << double(counter_nu_out_fv      * data_scale_factor  ) << std::endl;
    std::cout << " Number of Cosmic        : " << counter_cosmic         << "\t \t " << double(counter_cosmic         * data_scale_factor  ) << std::endl;
    std::cout << " Number of Numu CC       : " << counter_numu_cc        << "\t \t " << double(counter_numu_cc        * data_scale_factor  ) << std::endl;
    std::cout << " Number of Numu CC Pi0   : " << counter_numu_cc_pi0    << "\t \t " << double(counter_numu_cc_pi0    * data_scale_factor  ) << std::endl;
    std::cout << " Number of NC            : " << counter_nc             << "\t \t " << double(counter_nc             * data_scale_factor  ) << std::endl;
    std::cout << " Number of NC Pi0        : " << counter_nc_pi0         << "\t \t " << double(counter_nc_pi0         * data_scale_factor  ) << std::endl;
    std::cout << " Number of Unmatched     : " << counter_unmatched      << "\t \t " << double(counter_unmatched      * data_scale_factor  ) << std::endl;
    std::cout << " Number of InTime Cosmics: " << double(counter_ext * (intime_scale_factor / data_scale_factor))
              << "\t \t " << double(counter_ext * intime_scale_factor) << std::endl;
    std::cout << " Number of Dirt          : " << double(counter_dirt * dirt_scale_factor / data_scale_factor)
              << "\t \t " << double (counter_dirt * dirt_scale_factor)<< std::endl;
    std::cout << "--------- Neutrinos Selected in Truth ----------" << std::endl;
    std::cout << " Nue CC QE               : " << counter_nue_cc_qe   << std::endl;
    std::cout << " Nue CC Res              : " << counter_nue_cc_res  << std::endl;
    std::cout << " Nue CC DIS              : " << counter_nue_cc_dis  << std::endl;
    std::cout << " Nue CC COH              : " << counter_nue_cc_coh  << std::endl;
    std::cout << " Nue CC MEC              : " << counter_nue_cc_mec  << std::endl;
    std::cout << " Numu CC QE              : " << counter_numu_cc_qe  << std::endl;
    std::cout << " Numu CC Res             : " << counter_numu_cc_res << std::endl;
    std::cout << " Numu CC DIS             : " << counter_numu_cc_dis << std::endl;
    std::cout << " Numu CC COH             : " << counter_numu_cc_coh << std::endl;
    std::cout << " Numu CC MEC             : " << counter_numu_cc_mec << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " Tot Nue                 : " << counter_nue_cc_qe    + counter_nue_cc_res   + counter_nue_cc_dis   +  counter_nue_cc_coh   + counter_nue_cc_mec  << std::endl;
    std::cout << " Tot NuMu                : " << counter_numu_cc_qe   + counter_numu_cc_res  + counter_numu_cc_dis  +  counter_numu_cc_coh  + counter_numu_cc_mec << std::endl;
    std::cout << " Tot NC                  : " << counter_tot_nue_numu_nc << std::endl;
    std::cout << " Sum Neutrinos           : " <<   counter_nue_cc_qe   + counter_nue_cc_res  + counter_nue_cc_dis  +
                                                  + counter_nue_cc_coh  + counter_nue_cc_mec  + counter_numu_cc_qe  + counter_numu_cc_res +
                                                  + counter_numu_cc_dis + counter_numu_cc_coh + counter_numu_cc_mec + 
                                                  counter_tot_nue_numu_nc << std::endl;

    std::cout << "------------------------------------------------" << std::endl;
    const double efficiency = double(counter_nue_cc) / double(tot_true_infv_nues);
    const double purity     = double(counter_nue_cc) / double(counter);
    std::cout << " Efficiency       : " << "( " << counter_nue_cc << " / " << tot_true_infv_nues << " ) = " << efficiency << std::endl;
    std::cout << " Purity           : " << "( " << counter_nue_cc << " / " << counter           << " ) = " << purity << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " Total Nue Candidates in data : " << counter_data << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------