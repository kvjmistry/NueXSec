#include "../include/utility.h"

// -----------------------------------------------------------------------------
void utility::Initalise(){

    std::cout << "Initialising Utility Class..." << std::endl;

    std::string line;

    config_v.resize(k_config_MAX, 1.0);

    std::ifstream myfile ("config.txt");

    std::string varname;
    std::string value;
    
    if (myfile.is_open()) {

        // Loop over the config ist
        for (unsigned int p = 0; p < confignames.size(); p++){
            
            // Loop over lines in file
            while ( getline (myfile,line) ) {

                std::istringstream ss(line);
                ss >> varname >> value; 

                // Found the correct variation file 
                if (varname == confignames.at(p)) {
                    // std::cout << "Found match for: " << varname << " "<< std::stod(value) <<  std::endl;
                    config_v.at(p)= std::stod(value);
                    break;
                }
                
            }
        }

        myfile.close();
    }
    else std::cout << "Unable to open file, bad things are going to happen..." << std::endl; 

}
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
bool utility::GetHist(TFile* f, TH2D* &h, TString string){
    h = (TH2D*) f->Get(string);
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
void utility::Tabulate(bool inFV, std::string interaction, std::string classification, int type, std::vector<double> &counter_v, double weight) {

    if (type == k_mc){
        
        // Require in FV condition
        if (inFV){
            if (interaction == "nue_cc_qe")  counter_v.at(k_count_nue_cc_qe)  += weight; 
            if (interaction == "nue_cc_res") counter_v.at(k_count_nue_cc_res) += weight;
            if (interaction == "nue_cc_coh") counter_v.at(k_count_nue_cc_dis) += weight;
            if (interaction == "nue_cc_dis") counter_v.at(k_count_nue_cc_coh) += weight;
            if (interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_mec) += weight;

            if (interaction == "nue_bar_cc_qe")  counter_v.at(k_count_nuebar_cc_qe)  += weight; 
            if (interaction == "nue_bar_cc_res") counter_v.at(k_count_nuebar_cc_res) += weight;
            if (interaction == "nue_bar_cc_coh") counter_v.at(k_count_nuebar_cc_dis) += weight;
            if (interaction == "nue_bar_cc_dis") counter_v.at(k_count_nuebar_cc_coh) += weight;
            if (interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_mec) += weight;

            if (interaction == "numu_cc_qe"  )  counter_v.at(k_count_numu_cc_qe)  += weight;
            if (interaction == "numu_cc_res" )  counter_v.at(k_count_numu_cc_res) += weight;
            if (interaction == "numu_cc_coh" )  counter_v.at(k_count_numu_cc_dis) += weight;
            if (interaction == "numu_cc_dis" )  counter_v.at(k_count_numu_cc_coh) += weight;
            if (interaction == "numu_cc_mec" )  counter_v.at(k_count_numu_cc_mec) += weight;

            if (interaction == "numu_bar_cc_qe")  counter_v.at(k_count_numubar_cc_qe)  += weight;
            if (interaction == "numu_bar_cc_res") counter_v.at(k_count_numubar_cc_res) += weight;
            if (interaction == "numu_bar_cc_coh") counter_v.at(k_count_numubar_cc_dis) += weight;
            if (interaction == "numu_bar_cc_dis") counter_v.at(k_count_numubar_cc_coh) += weight;
            if (interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_mec) += weight;

            // These are all the nus, but now in the fv
            if (interaction == "nue_cc_qe" || interaction == "nue_cc_res" || interaction == "nue_cc_coh" || interaction == "nue_cc_dis" || interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_infv) += weight;
            if (interaction == "nue_bar_cc_qe" || interaction == "nue_bar_cc_res" || interaction == "nue_bar_cc_coh" || interaction == "nue_bar_cc_dis" || interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_infv) += weight;
            
            if (interaction == "numu_cc_qe" || interaction == "numu_cc_res" || interaction == "numu_cc_coh" || interaction == "numu_cc_dis" || interaction == "numu_cc_mec") counter_v.at(k_count_numu_cc_infv) += weight;
            if (interaction == "numu_bar_cc_qe" || interaction == "numu_bar_cc_res" || interaction == "numu_bar_cc_coh" || interaction == "numu_bar_cc_dis" || interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_infv) += weight;
        }

        // These are all the nus, but now in the cryostat volume
        if (interaction == "nue_cc_qe" || interaction == "nue_cc_res" || interaction == "nue_cc_coh" || interaction == "nue_cc_dis" || interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_incryo) += weight;
        if (interaction == "nue_bar_cc_qe" || interaction == "nue_bar_cc_res" || interaction == "nue_bar_cc_coh" || interaction == "nue_bar_cc_dis" || interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_incryo) += weight;
        
        if (interaction == "numu_cc_qe" || interaction == "numu_cc_res" || interaction == "numu_cc_coh" || interaction == "numu_cc_dis" || interaction == "numu_cc_mec") counter_v.at(k_count_numu_cc_incryo) += weight;
        if (interaction == "numu_bar_cc_qe" || interaction == "numu_bar_cc_res" || interaction == "numu_bar_cc_coh" || interaction == "numu_bar_cc_dis" || interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_incryo) += weight;
        
        // Classification
        if (classification == "nue_cc")       counter_v.at(k_count_nue_cc)       += weight;
        if (classification == "nue_cc_mixed") counter_v.at(k_count_nue_cc_mixed) += weight;
        if (classification == "nu_out_fv")    counter_v.at(k_count_nu_out_fv)    += weight;
        if (classification == "nc")           counter_v.at(k_count_nc)           += weight;
        if (classification == "nc_pi0")       counter_v.at(k_count_nc_pi0)       += weight;
        if (classification == "numu_cc")      counter_v.at(k_count_numu_cc)      += weight;
        if (classification == "numu_cc_pi0")  counter_v.at(k_count_numu_cc_pi0)  += weight;
        if (classification == "cosmic")       counter_v.at(k_count_cosmic)       += weight;
        if (classification == "unmatched")    counter_v.at(k_count_unmatched)    += weight;

        // Total selected MC events
        counter_v.at(k_count_total_mc) += weight;
    
    }
    else if (type == k_data) {
        counter_v.at(k_count_data) += weight; // The weight **should** always be 1 for these
    }
    else if (type == k_ext){
        counter_v.at(k_count_ext)  += weight; // The weight **should** always be 1 for these -- maybe not this if we decide to use the BNB stream
    }
    else if (type == k_dirt){
        counter_v.at(k_count_dirt) += weight;
    }
    else {

        std::cout << "unkown type specified!!!  " << __PRETTY_FUNCTION__ << std::endl;
    }
    
}
// -----------------------------------------------------------------------------
void utility::PrintInfo(std::vector<double> c_v, double intime_scale_factor, double mc_scale_factor, double dirt_scale_factor, std::string cut_name, double tot_true_infv_nues, double &efficiency, double &purity) {

    // c_v is short for counter vec here!

    // Sum of selected mc, dirt and ext. The dirt and ext are scaled to the MC POT
    double sum_mc_dirt_ext = c_v.at(k_count_total_mc)+ (c_v.at(k_count_ext) * (intime_scale_factor / mc_scale_factor)) + (c_v.at(k_count_dirt) * (dirt_scale_factor / mc_scale_factor));

    std::cout << "\n------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "\n\033[0;33m <" << cut_name << "> \033[0m" << std::endl;
    std::cout << "                    Scaled to MC POT | Scaled to Data POT | Unscaled" << std::endl;
    std::cout << " Total Candidate Nue     : " << sum_mc_dirt_ext                 << "\t "    << double(sum_mc_dirt_ext                * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Nue CC        : " << c_v.at(k_count_nue_cc)          << "\t \t " << double(c_v.at(k_count_nue_cc)         * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Nue CC Mixed  : " << c_v.at(k_count_nue_cc_mixed)    << "\t \t " << double(c_v.at(k_count_nue_cc_mixed)   * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Nu out FV     : " << c_v.at(k_count_nu_out_fv)       << "\t \t " << double(c_v.at(k_count_nu_out_fv)      * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Cosmic        : " << c_v.at(k_count_cosmic)          << "\t \t " << double(c_v.at(k_count_cosmic)         * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Numu CC       : " << c_v.at(k_count_numu_cc)         << "\t \t " << double(c_v.at(k_count_numu_cc)        * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Numu CC Pi0   : " << c_v.at(k_count_numu_cc_pi0)     << "\t \t " << double(c_v.at(k_count_numu_cc_pi0)    * mc_scale_factor  ) << std::endl;
    std::cout << " Number of NC            : " << c_v.at(k_count_nc)              << "\t \t " << double(c_v.at(k_count_nc)             * mc_scale_factor  ) << std::endl;
    std::cout << " Number of NC Pi0        : " << c_v.at(k_count_nc_pi0)          << "\t \t " << double(c_v.at(k_count_nc_pi0)         * mc_scale_factor  ) << std::endl;
    std::cout << " Number of Unmatched     : " << c_v.at(k_count_unmatched)       << "\t \t " << double(c_v.at(k_count_unmatched)      * mc_scale_factor  ) << std::endl;
    
    std::cout << " Number of InTime Cosmics: " << double(c_v.at(k_count_ext) * (intime_scale_factor / mc_scale_factor))
              << "\t " << double(c_v.at(k_count_ext) * intime_scale_factor) << "\t " << c_v.at(k_count_ext) << std::endl;
    
    std::cout << " Number of Dirt          : " << double(c_v.at(k_count_dirt) * dirt_scale_factor / mc_scale_factor)
              << "\t "                         << double(c_v.at(k_count_dirt) * dirt_scale_factor) << "\t " << c_v.at(k_count_dirt) << std::endl;

    std::cout << "----------- Neutrinos in FV Truth -------------" << std::endl;
    std::cout << " Nue CC QE    : " << c_v.at(k_count_nue_cc_qe)   <<  "   Nuebar CC QE    : " << c_v.at(k_count_nuebar_cc_qe)  << std::endl;
    std::cout << " Nue CC Res   : " << c_v.at(k_count_nue_cc_res)  <<  "   Nuebar CC Res   : " << c_v.at(k_count_nuebar_cc_res) << std::endl;
    std::cout << " Nue CC DIS   : " << c_v.at(k_count_nue_cc_dis)  <<  "   Nuebar CC DIS   : " << c_v.at(k_count_nuebar_cc_dis) << std::endl;
    std::cout << " Nue CC COH   : " << c_v.at(k_count_nue_cc_coh)  <<  "   Nuebar CC COH   : " << c_v.at(k_count_nuebar_cc_coh) << std::endl;
    std::cout << " Nue CC MEC   : " << c_v.at(k_count_nue_cc_mec)  <<  "   Nuebar CC MEC   : " << c_v.at(k_count_nuebar_cc_mec) << std::endl;
    std::cout << std::endl;
    std::cout << " Numu CC QE   : " << c_v.at(k_count_numu_cc_qe)  <<  "   Numubar CC QE   : " << c_v.at(k_count_numubar_cc_qe) << std::endl;
    std::cout << " Numu CC Res  : " << c_v.at(k_count_numu_cc_res) <<  "   Numubar CC Res  : " << c_v.at(k_count_numubar_cc_res)<< std::endl;
    std::cout << " Numu CC DIS  : " << c_v.at(k_count_numu_cc_dis) <<  "   Numubar CC DIS  : " << c_v.at(k_count_numubar_cc_dis)<< std::endl;
    std::cout << " Numu CC COH  : " << c_v.at(k_count_numu_cc_coh) <<  "   Numubar CC COH  : " << c_v.at(k_count_numubar_cc_coh)<< std::endl;
    std::cout << " Numu CC MEC  : " << c_v.at(k_count_numu_cc_mec) <<  "   Numubar CC MEC  : " << c_v.at(k_count_numubar_cc_mec)<< std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " Tot Nue in FV                 : " << c_v.at(k_count_nue_cc_infv)     << "   Tot Nue in Cryo                 : " << c_v.at(k_count_nue_cc_incryo)     << std::endl;
    std::cout << " Tot Nuebar in FV              : " << c_v.at(k_count_nuebar_cc_infv)  << "   Tot Nuebar in Cryo              : " << c_v.at(k_count_nuebar_cc_incryo) << std::endl;
    std::cout << " Tot NuMu in FV                : " << c_v.at(k_count_numu_cc_infv)    << "   Tot NuMu in Cryo                : " << c_v.at(k_count_numu_cc_incryo)<< std::endl;
    std::cout << " Tot NuMubar in FV             : " << c_v.at(k_count_numubar_cc_infv) << "   Tot NuMubar in Cryo             : " << c_v.at(k_count_numubar_cc_incryo)<< std::endl;
    // std::cout << " Tot NC                  : " << counter_tot_nue_numu_nc << std::endl;
    // std::cout << " Sum Neutrinos           : " <<  << std::endl;

    std::cout << "------------------------------------------------" << std::endl;
    efficiency = double(c_v.at(k_count_nue_cc)) / double(tot_true_infv_nues);
    purity     = double(c_v.at(k_count_nue_cc)) / double(sum_mc_dirt_ext);
    std::cout << " Efficiency       : " << "( " << c_v.at(k_count_nue_cc) << " / " << tot_true_infv_nues << " ) = " << efficiency << std::endl;
    std::cout << " Purity           : " << "( " << c_v.at(k_count_nue_cc) << " / " << sum_mc_dirt_ext           << " ) = " << purity << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << " Total Nue Candidates in data : " << c_v.at(k_count_data) << std::endl;
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

    // From beam to detector rotation matrix
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam
    // RotDet2Beam.Invert(); // Invert back to the beam to det

    // Rotate to beam coords
    BeamCoords = RotDet2Beam * detxyz;

    TVector3 beam_dir = {0 , 0 , 1};
    double theta = BeamCoords.Angle(beam_dir) * 180 / 3.1415926;

    // std::cout << theta << std::endl;

    return theta;
}
// -----------------------------------------------------------------------------
