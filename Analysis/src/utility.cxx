#include "../include/utility.h"


// -----------------------------------------------------------------------------
std::vector<double> utility::configure(
        double _Run1_MC_POT,
        double _Run1_Dirt_POT,
        double _Run1_Data_POT,
        double _Run1_Data_trig,
        double _Run1_EXT_trig,
        double _x1,
        double _x2,
        double _y1,
        double _y2,
        double _z1,
        double _z2
        ) {
    std::vector<double> config;
    config.resize(k_config_MAX, 0);

    config.at(k_config_Run1_MC_POT)     = _Run1_MC_POT;
    config.at(k_config_Run1_Dirt_POT)   = _Run1_Dirt_POT;
    config.at(k_config_Run1_Data_POT)   = _Run1_Data_POT;
    config.at(k_config_Run1_Data_trig)  = _Run1_Data_trig;
    config.at(k_config_Run1_EXT_trig)   = _Run1_EXT_trig;

    config.at(k_config_x1)  = _x1;
    config.at(k_config_x2)  = _x2;
    config.at(k_config_y1)  = _y1;
    config.at(k_config_y2)  = _y2;
    config.at(k_config_z1)  = _z1;
    config.at(k_config_z2)  = _z2;
    

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
        if (verbose) std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
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
void utility::Tabulate(std::string interaction, std::string classification, int type, std::vector<double> &counter_v, double weight) {

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
        double nue_cc        = 0.0;
        double nue_cc_mixed  = 0.0;
        double cosmic        = 0.0;
        double numu_cc       = 0.0;
        double numu_cc_pi0   = 0.0;
        double nu_out_fv     = 0.0;
        double nc            = 0.0;
        double nc_pi0        = 0.0;
        double unmatched     = 0.0;
        double total         = 0.0;

        // Total nue or nuebar interaction types
        int only_nue_cc     = 0;
        int only_nue_bar_cc = 0;

        if (interaction == "nue_cc_qe")  nue_cc_qe  += weight; 
        if (interaction == "nue_cc_res") nue_cc_res += weight;
        if (interaction == "nue_cc_coh") nue_cc_coh += weight;
        if (interaction == "nue_cc_dis") nue_cc_dis += weight;
        if (interaction == "nue_cc_mec") nue_cc_mec += weight;

        if (interaction == "nue_bar_cc_qe")  nue_bar_cc_qe  += weight; 
        if (interaction == "nue_bar_cc_res") nue_bar_cc_res += weight;
        if (interaction == "nue_bar_cc_coh") nue_bar_cc_coh += weight;
        if (interaction == "nue_bar_cc_dis") nue_bar_cc_dis += weight;
        if (interaction == "nue_bar_cc_mec") nue_bar_cc_mec += weight;

        if (interaction == "numu_cc_qe"  || interaction == "numu_bar_cc_qe")  numu_cc_qe  += weight;
        if (interaction == "numu_cc_res" || interaction == "numu_bar_cc_res") numu_cc_res += weight;
        if (interaction == "numu_cc_coh" || interaction == "numu_bar_cc_coh") numu_cc_coh += weight;
        if (interaction == "numu_cc_dis" || interaction == "numu_bar_cc_dis") numu_cc_dis += weight;
        if (interaction == "numu_cc_mec" || interaction == "numu_bar_cc_mec") numu_cc_mec += weight;
        
        // Check if string contains "NC"
        if (interaction.find("nc") != std::string::npos) {
            tot_nue_numu_nc += weight;
        }

        if (classification == "nue_cc")       nue_cc       += weight;
        if (classification == "nue_cc_mixed") nue_cc_mixed += weight;
        if (classification == "nu_out_fv")    nu_out_fv    += weight;
        if (classification == "nc")           nc           += weight;
        if (classification == "nc_pi0")       nc_pi0       += weight;
        if (classification == "numu_cc")      numu_cc      += weight;
        if (classification == "numu_cc_pi0")  numu_cc_pi0  += weight;
        if (classification == "cosmic")       cosmic       += weight;
        if (classification == "unmatched")    unmatched    += weight;

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
        counter_v.at(k_count_total_nue_cc_qe)  += total_nue_cc_qe ;
        counter_v.at(k_count_total_nue_cc_res) += total_nue_cc_res;
        counter_v.at(k_count_total_nue_cc_dis) += total_nue_cc_dis;
        counter_v.at(k_count_total_nue_cc_coh) += total_nue_cc_coh;
        counter_v.at(k_count_total_nue_cc_mec) += total_nue_cc_mec;

        counter_v.at(k_count_only_nue_cc)      += only_nue_cc    ;
        counter_v.at(k_count_only_nue_bar_cc)  += only_nue_bar_cc;

        counter_v.at(k_count_numu_cc_qe)       += numu_cc_qe ;
        counter_v.at(k_count_numu_cc_res)      += numu_cc_res;
        counter_v.at(k_count_numu_cc_dis)      += numu_cc_dis;
        counter_v.at(k_count_numu_cc_coh)      += numu_cc_coh;
        counter_v.at(k_count_numu_cc_mec)      += numu_cc_mec;

        counter_v.at(k_count_tot_nue_numu_nc)  += tot_nue_numu_nc;
        
        counter_v.at(k_count_nue_cc)           += nue_cc      ;
        counter_v.at(k_count_nue_cc_mixed)     += nue_cc_mixed;
        counter_v.at(k_count_nu_out_fv)        += nu_out_fv   ;
        counter_v.at(k_count_cosmic)           += cosmic      ;
        counter_v.at(k_count_numu_cc)          += numu_cc     ;
        counter_v.at(k_count_numu_cc_pi0)      += numu_cc_pi0 ;
        counter_v.at(k_count_nc)               += nc          ;
        counter_v.at(k_count_nc_pi0)           += nc_pi0      ;
        counter_v.at(k_count_unmatched)        += unmatched   ;
        
        counter_v.at(k_count_total)            += total;
    
    }
    else if (type == k_data) {
        counter_v.at(k_count_data) += 1;
    }
    else if (type == k_ext){
        counter_v.at(k_count_ext)  += 1;
    }
    else if (type == k_dirt){
        counter_v.at(k_count_dirt) += weight;
    }
    else {

        std::cout << "unkown type specified!!!  " << __PRETTY_FUNCTION__ << std::endl;
    }
    
}
// -----------------------------------------------------------------------------
void utility::PrintInfo(std::vector<double> counter_v, double intime_scale_factor, double mc_scale_factor, double dirt_scale_factor, std::string cut_name, double tot_true_infv_nues, double &efficiency, double &purity) {

    double counter_nue_cc_qe      = counter_v.at(k_count_total_nue_cc_qe);
    double counter_nue_cc_res     = counter_v.at(k_count_total_nue_cc_res);
    double counter_nue_cc_dis     = counter_v.at(k_count_total_nue_cc_dis);
    double counter_nue_cc_coh     = counter_v.at(k_count_total_nue_cc_coh);
    double counter_nue_cc_mec     = counter_v.at(k_count_total_nue_cc_mec);
    
    double counter_numu_cc_qe     = counter_v.at(k_count_numu_cc_qe);
    double counter_numu_cc_res    = counter_v.at(k_count_numu_cc_res);
    double counter_numu_cc_dis    = counter_v.at(k_count_numu_cc_dis);
    double counter_numu_cc_coh    = counter_v.at(k_count_numu_cc_coh);
    double counter_numu_cc_mec    = counter_v.at(k_count_numu_cc_mec);

    double counter_tot_nue_numu_nc = counter_v.at(k_count_tot_nue_numu_nc);

    double counter_nue_cc         = counter_v.at(k_count_nue_cc);
    double counter_nue_cc_mixed   = counter_v.at(k_count_nue_cc_mixed);
    double counter_nu_out_fv      = counter_v.at(k_count_nu_out_fv);
    double counter_cosmic         = counter_v.at(k_count_cosmic);
    double counter_numu_cc        = counter_v.at(k_count_numu_cc);
    double counter_numu_cc_pi0    = counter_v.at(k_count_numu_cc_pi0);
    double counter_nc             = counter_v.at(k_count_nc);
    double counter_nc_pi0         = counter_v.at(k_count_nc_pi0);
    double counter_unmatched      = counter_v.at(k_count_unmatched);
    double counter                = counter_v.at(k_count_total);
    
    double counter_data           = counter_v.at(k_count_data);
    double counter_ext            = counter_v.at(k_count_ext);
    double counter_dirt           = counter_v.at(k_count_dirt);

    counter = counter + (counter_ext * (intime_scale_factor / mc_scale_factor)) + (counter_dirt * (dirt_scale_factor / mc_scale_factor)); // Add in the EXT and dirt counters with the scalings to MC POT

    std::cout << "\n------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "\n\033[0;33m <" << cut_name << "> \033[0m" << std::endl;
    std::cout << "                    Scaled to MC POT | Scaled to Data POT | Unscaled" << std::endl;
    std::cout << " Total Candidate Nue     : " << counter                << "\t " << double(counter                   * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Nue CC        : " << counter_nue_cc         << "\t \t " << double(counter_nue_cc         * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Nue CC Mixed  : " << counter_nue_cc_mixed   << "\t \t " << double(counter_nue_cc_mixed   * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Nu out FV     : " << counter_nu_out_fv      << "\t \t " << double(counter_nu_out_fv      * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Cosmic        : " << counter_cosmic         << "\t \t " << double(counter_cosmic         * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Numu CC       : " << counter_numu_cc        << "\t \t " << double(counter_numu_cc        * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Numu CC Pi0   : " << counter_numu_cc_pi0    << "\t \t " << double(counter_numu_cc_pi0    * mc_scale_factor  ) << std::endl;
    std::cout << " Number of NC            : " << counter_nc             << "\t \t " << double(counter_nc             * mc_scale_factor  ) << std::endl;
    std::cout << " Number of NC Pi0        : " << counter_nc_pi0         << "\t \t " << double(counter_nc_pi0         * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Unmatched     : " << counter_unmatched      << "\t \t " << double(counter_unmatched      * mc_scale_factor  ) << std::endl;
    
    std::cout << " Number of InTime Cosmics: " << double(counter_ext * (intime_scale_factor / mc_scale_factor))
              << "\t " << double(counter_ext * intime_scale_factor) << "\t " << counter_ext << std::endl;
    
    std::cout << " Number of Dirt          : " << double(counter_dirt * dirt_scale_factor / mc_scale_factor)
              << "\t "                         << double(counter_dirt * dirt_scale_factor) << "\t " << counter_dirt << std::endl;

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
    efficiency = double(counter_nue_cc) / double(tot_true_infv_nues);
    purity     = double(counter_nue_cc) / double(counter);
    std::cout << " Efficiency       : " << "( " << counter_nue_cc << " / " << tot_true_infv_nues << " ) = " << efficiency << std::endl;
    std::cout << " Purity           : " << "( " << counter_nue_cc << " / " << counter           << " ) = " << purity << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " Total Nue Candidates in data : " << counter_data << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------
double utility::GetTheta(double px, double py, double pz){

    // Variables
    TRotation RotDet2Beam;             // Rotations
    TVector3  detxyz, BeamCoords;      // Translations
    std::vector<double> rotmatrix;     // Inputs

    // input detector coordinates to translate
    detxyz = {px, py, pz};     

    // From detector to beam coords
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now beam to det
    // RotDet2Beam.Invert(); // Invert back to the det to beam rot matrix

    TVector3 Trans_Targ2Det_beam = { 5502, 7259, 67270}; //cm -- in beam coords

    // Rotate to beam coords
    BeamCoords = RotDet2Beam * detxyz + Trans_Targ2Det_beam;

    TVector3 beam_dir = {0 , 0 , 1};
    double theta = BeamCoords.Angle(beam_dir) * 180 / 3.1415926;

    // std::cout << theta << std::endl;

    return theta;
}
// -----------------------------------------------------------------------------
