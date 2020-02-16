#include "../include/utility.h"


// -------------------------------------------------------------------------
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
// -------------------------------------------------------------------------
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
// -------------------------------------------------------------------------
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
// -------------------------------------------------------------------------
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
// -------------------------------------------------------------------------
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
// -------------------------------------------------------------------------

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
