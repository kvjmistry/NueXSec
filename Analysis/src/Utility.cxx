#include "../include/Utility.h"

// -----------------------------------------------------------------------------
void Utility::Initalise(int argc, char *argv[], std::string usage,std::string usage2, std::string usage3 ){

    std::cout << "Initialising Utility Class..." << std::endl;

    // Loop over input arguments
    for (int i =1; i < argc; i++){
        auto const arg = argv[i];
        //std::cout << arg << std::endl; // This is for debugging
        
        // Slim input
        if (strcmp(arg, "--slim") == 0) {
            std::cout << "Running with slim mode"<< std::endl;
            slim = true;
            std::cout << red << " *** \t Running with Slimmed Selection (no histograms will be made)\t *** " << reset << std::endl;
        }

        // Histogram Mode
        if (strcmp(arg, "--hist") == 0) {
            std::cout << yellow << "Making Histograms, file to make histograms with: "<< argv[i+1] << yellow <<std::endl;
            make_histos = true;
            run_selection = false; // switch this bool out
            hist_file_name = argv[i+1];
        }

        // XSec mode
        if (strcmp(arg, "--xsec") == 0) {
            std::cout << yellow << "Calculating the Cross Section with input file name: " << argv[i+1] << yellow << std::endl;
            calc_cross_sec = true;
            run_selection = false; // switch this bool out
            tree_file_name = argv[i+1];
        }
        
        // MC file
        if (strcmp(arg, "--mc") == 0) {
            // std::cout << "Running with MC file: " << argv[i+1] << std::endl;
            mc_file_name = argv[i+1];
        }

        // Overwrite output mc file name
        if (strcmp(arg, "--mc_out") == 0) {
            std::cout << "New Output MC File name: " << argv[i+1] << std::endl;
            mc_file_name_out = argv[i+1];
        }

        // EXT file
        if (strcmp(arg, "--ext") == 0){
            // std::cout << "Running with EXT file: " << argv[i+1] << std::endl;
            ext_file_name = argv[i+1];
        }

        // Overwrite output ext file name
        if (strcmp(arg, "--ext_out") == 0) {
            std::cout << "New Output EXT File name: " << argv[i+1] << std::endl;
            ext_file_name_out = argv[i+1];
        }

        // Data file
        if (strcmp(arg, "--data") == 0){
            // std::cout << "Running with Data file: " << argv[i+1] << std::endl;
            data_file_name = argv[i+1];
        }
        
        // Overwrite output data file name
        if (strcmp(arg, "--data_out") == 0) {
            std::cout << "New Output Data File name: " << argv[i+1] << std::endl;
            data_file_name_out = argv[i+1];
        }

        // Dirt file overlay
        if (strcmp(arg, "--dirt") == 0){
            // std::cout << "Running with Dirt file: " << argv[i+1] << std::endl;
            dirt_file_name = argv[i+1];
        }

        // Overwrite output dirt file name
        if (strcmp(arg, "--dirt_out") == 0) {
            std::cout << "New Output Dirt File name: " << argv[i+1] << std::endl;
            dirt_file_name_out = argv[i+1];
        }

        // GENIE Tune Weight Settings
        if (strcmp(arg, "--weight_tune") == 0){
            std::cout << "Running with GENIE Tune mode: " << argv[i+1] << std::endl;
            _weight_tune = atoi(argv[i+1]);
        }
        
        // PPFX CV Weight Settings
        if (strcmp(arg, "--weight_ppfx") == 0){
            std::cout << "Running with PPFX CV mode: " << argv[i+1] << std::endl;
            _weight_ppfx = atoi(argv[i+1]);
        }

        // Dirt Weight Settings
        if (strcmp(arg, "--weight_dirt") == 0){
            std::cout << "Running with Dirt mode: " << argv[i+1] << std::endl;
            _weight_dirt = atoi(argv[i+1]);
        }

        // EXT Weight Settings
        if (strcmp(arg, "--weight_ext") == 0){
            std::cout << "Running with EXT mode: " << argv[i+1] << std::endl;
            _weight_ext = atoi(argv[i+1]);
        }

        // pi0 Weight Settings
        if (strcmp(arg, "--weight_pi0") == 0){
            std::cout << "Running with pi0 mode: " << argv[i+1] << std::endl;
            _pi0_correction = atoi(argv[i+1]);
        }

        // Whats the verbose?
        if (strcmp(arg, "-v") == 0 || strcmp(arg, "--verbose") == 0 || strcmp(arg, "--v") == 0){
            std::cout << "Setting Verbose Level to : " << argv[i+1] << std::endl;
            verbose = atoi(argv[i+1]);
        }

        // Max number of events specified?
        if (strcmp(arg, "-n") == 0 || strcmp(arg, "--n") == 0){
            std::cout << red << "Running with a maximum of : " << argv[i+1] << " events" << reset <<std::endl;
            num_events = atoi(argv[i+1]);
        }

        // Set the run period
        if (strcmp(arg, "--run") == 0){
            // std::cout << "Setting the run period as : run" << argv[i+1] <<std::endl;
            run_period = argv[i+1];
        }

        // Help!
        if (strcmp(arg, "--h") == 0 || strcmp(arg, "-h") == 0|| strcmp(arg, "--help") == 0 || strcmp(arg, "--usage") == 0){
            std::cout << usage <<  usage2 << usage3 << std::endl; 
            exit(1);
        }

        // Area Normalise the histograms
        if (strcmp(arg, "--area") == 0) {
            std::cout << "Area Normalising the histograms"<< std::endl;
            area_norm = true;
        }

        // Variation file
        if (strcmp(arg, "--var") == 0){
            std::cout << magenta << "Using Variation: " << argv[i+2] << reset << std::endl;
            
            variation = argv[i+2];
            mc_file_name = argv[i+1];
            
            mc_file_name_out      = Form("nuexsec_mc_run%s_%s.root", run_period, variation);
            mc_tree_file_name_out = Form("nuexsec_selected_tree_mc_run%s_%s.root", run_period, variation);
            
            std::cout  <<yellow << "Output filename will be overwritten with name: "      << mc_file_name_out << reset <<  std::endl;
            std::cout <<yellow << "Output tree filename will be overwritten with name: " << mc_tree_file_name_out << reset <<  std::endl;
            
            overwritePOT = true;
        }

        // Systematics
        if (strcmp(arg, "--sys") == 0){
            std::cout << "Using Systematics plotting code with mode: " << argv[i+1] << std::endl;
            run_sys = true;
            run_selection = false;
            sysmode = argv[i+1];
        }

        // Plot systematics on the CV histograms
        if (strcmp(arg, "--plotsys") == 0){
            std::cout << "Plotting systematic uncertainties on the CV Histograms for: " << argv[i+1] << std::endl;
            plot_sys_uncertainty = true;
            sysplot = argv[i+1];
        }

        // Cross-Section
        if (strcmp(arg, "--xsecmode") == 0){
            std::cout << "Using Cross-Section code with mode: " << argv[i+1] << std::endl;
            xsecmode = argv[i+1];
            xsec_labels = argv[i+2];
            xsec_rw_mode = argv[i+3];
        }

        // Utility Plotter
        if (strcmp(arg, "--uplot") == 0){
            std::cout << "Using Utility plotting code with mode: " << argv[i+1] << std::endl;
            run_uplot = true;
            run_selection = false;
            uplotmode = argv[i+1];
        }

        // Only run the print function
        if (strcmp(arg, "--printonly") == 0) {
            run_selection = false;
            print = true;
        }

        // Print all
        if (strcmp(arg, "--printall") == 0) {
            print = true;
            print_mc                 = true;
            print_data               = true;
            print_ext                = true;
            print_dirt               = true;
        }
        
        // Print MC
        if (strcmp(arg, "--printmc") == 0) {
            print = true;
            print_mc = true;
        }
        
        // Print Data
        if (strcmp(arg, "--printdata") == 0) {
            print = true;
            print_data = true;
        }

        // Print ext
        if (strcmp(arg, "--printext") == 0) {
            print = true;
            print_ext = true;
        }

        // Print Dirt
        if (strcmp(arg, "--printdirt") == 0) {
            print = true;
            print_dirt = true;
        }

        // Intrinsic nue Mode
        if (strcmp(arg, "--intrinsic") == 0) {
            std::cout << magenta << "Using Selection with intrinsic nue setting: " << argv[i+1] << reset <<std::endl;
            intrinsic_mode = argv[i+1];
        }

        // Use the gpvm file paths
        if (strcmp(arg, "--gpvm") == 0) {
            use_gpvm = true;
        }
   
    }

    // Finished getting all the input parameters to the code
    // ----------------------------------------------------------------------------

    std::string variation_str = variation;

    // Check the variaition mode setting and see if we want to overwrite the POT
    if (overwritePOT){

        // Loop over the POT config names and overwrite the name of the CV MC POT
        for (unsigned int p = 0; p < confignames.size(); p++){
            
            std::string match_name = Form("Run%s_MC_POT", run_period);
            // std::string match_name = "Run1_MC_POT";
            
            // If matched then overwrite the POT config for the MC to the variation
            if (confignames.at(p) == match_name){
                confignames[p] = match_name + "_" + variation_str;
                std::cout << red  <<"New MC POT config to search for is: " << confignames.at(p) << reset << std::endl;
            }
        }

        // We also want to set the intrinsic nue POT so the weighting is done correctly
        // Loop over the POT config names and overwrite the name of the CV MC POT
        for (unsigned int p = 0; p < confignames.size(); p++){
            
            std::string match_name = Form("Run%s_Intrinsic_POT", run_period);
            // std::string match_name = "Run1_Intrinsic_POT";
            
            // If matched then overwrite the POT config for the MC to the variation
            if (confignames.at(p) == match_name){
                confignames[p] = match_name + "_" + variation_str;
                std::cout << red << "New MC POT config to search for is: " << confignames.at(p) << reset <<std::endl;
            }
        }
    
    
    }

    // Get the congigureation parameters
    std::string line;

    config_v.resize(k_config_MAX, 1.0);

    std::string varname;
    std::string value;
    
    // Loop over the config ist
    for (unsigned int p = 0; p < confignames.size(); p++){

        std::ifstream myfile ("config.txt");

        if (myfile.is_open()) {
            
            // std::cout << confignames.at(p) <<  std::endl;

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
        else std::cout << "Unable to open file, bad things are going to happen..." << std::endl; 

        myfile.close();
    }

    // check if the POT from the input file and in config.txt match
    // CheckPOT();

    // Now set the weight configurations
    
    // GENIE Tune
    if (_weight_tune == 0){
        weight_tune = false;
        std::cout << "GENIE Tune is turned off" << std::endl;
    }
    else {
        weight_tune = true;
        std::cout << "Using GENIE Tune"<< std::endl;
    }

    // PPFX CV
    if (_weight_ppfx == 0){
        weight_ppfx = false;
        std::cout << "PPFX CV Correction is turned off" << std::endl;
    }
    else {
        weight_ppfx = true;
        std::cout << "Using PPFX CV Correction"<< std::endl;
    }

    // Dirt correction
    if (_weight_dirt == 0){
        weight_dirt = false;
        std::cout << "Dirt Correction is turned off" << std::endl;
    }
    else {
        weight_dirt = true;
        std::cout << "Using Correction factor for dirt"<< std::endl;
    }

    // EXT correction
    if (_weight_ext == 0){
        weight_ext = false;
        std::cout << "EXT Correction is turned off" << std::endl;
    }
    else {
        weight_ext = true;
        std::cout << "Using Correction factor for ext"<< std::endl;
    }

    // Pi0 correction
    pi0_correction = _pi0_correction;
    if (pi0_correction == 0){
        std::cout << blue << "pi0 Correction is turned off" << reset << std::endl;
    }
    else if (pi0_correction == 1){
         std::cout << "Using "<< blue << "normalisation factor" << reset<< " to correct pi0" << std::endl;
    }
    else {
        std::cout << "Using :"<< blue << "energy dependent scaling" << reset << "to correct pi0" << std::endl;
    }

    // Set the scale factors
    if (strcmp(run_period, "1") == 0){
        mc_scale_factor     = config_v.at(k_Run1_Data_POT)  / config_v.at(k_Run1_MC_POT);
        dirt_scale_factor   = 0.45*config_v.at(k_Run1_Data_POT)  / config_v.at(k_Run1_Dirt_POT);
        ext_scale_factor    = config_v.at(k_Run1_Data_trig) / config_v.at(k_Run1_EXT_trig);
        intrinsic_weight    = config_v.at(k_Run1_MC_POT)    / config_v.at(k_Run1_Intrinsic_POT);
    }
    else if (strcmp(run_period, "3") == 0){
        mc_scale_factor     = config_v.at(k_Run3_Data_POT)  / config_v.at(k_Run3_MC_POT);
        dirt_scale_factor   = 0.45*config_v.at(k_Run3_Data_POT)  / config_v.at(k_Run3_Dirt_POT);
        ext_scale_factor    = 0.95*config_v.at(k_Run3_Data_trig) / config_v.at(k_Run3_EXT_trig);
    }
    else {
        std::cout << "Error Krish... You havent specified the run period!" << std::endl;
        std::cout << usage <<  usage2 << usage3 << std::endl;
        exit(1);
    }

    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << red << "Scale Factors:\n" << green <<
    "MC Scale factor:   "   << mc_scale_factor     << "\n" <<
    "Dirt Scale factor: "   << dirt_scale_factor   << "\n" <<
    "EXT Scale factor:  "   << ext_scale_factor << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    // Display the intrinsic nue scale factor
    if (std::string(intrinsic_mode) == "intrinsic"){
        std::cout << blue << "Intrinsic Factor:\n" << intrinsic_weight << std::endl;
        std::cout << "-------------------------------" << reset << std::endl;
    }

    if (use_gpvm)
        std::cout << red << "Using gpvm environment, all paths will be set compatible with running on a gpvm"<< reset << std::endl;
    else 
        std::cout << red << "Using local environment, all paths will be set compatible with running on Krish's computer"<< reset << std::endl;
    
    
}
// -----------------------------------------------------------------------------
bool Utility::GetFile(TFile* &f, TString string){
    
    if (string == "empty") return false;
    
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
bool Utility::GetTree(TFile* f, TTree* &T, TString string){
    f->cd();
    T = (TTree*) f->Get(string);
    if (T == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis tree might not exist in the file\n" << std::endl;
        // exit(1);
        return true;
    }
    else {
        return false;
    }
}
// -----------------------------------------------------------------------------
bool Utility::GetHist(TFile* f, TH1D* &h, TString string){
    f->cd();
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
bool Utility::GetHist(TFile* f, TH2D* &h, TString string){
    f->cd();
    h = (TH2D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
void Utility::CheckWeight(double &weight){

    // Infinate weight
    if (std::isinf(weight))           weight = 1.0;
    
    // NAN weight 
    else if (std::isnan(weight) == 1) weight = 1.0;
    
    // Large weight
    else if (weight > 30)             weight = 1.0;
    
    // Negative weight
    else if (weight < 0)              weight = 1.0;
    
    else {
        // Do nothing to the weight
    }

}
// -----------------------------------------------------------------------------
void Utility::CheckWeight(float &weight){

    // Infinate weight
    if (std::isinf(weight))           weight = 1.0;
    
    // NAN weight 
    else if (std::isnan(weight) == 1) weight = 1.0;
    
    // Large weight
    else if (weight > 30)             weight = 1.0;
    
    // Negative weight
    else if (weight < 0)              weight = 1.0;
    
    else {
        // Do nothing to the weight
    }

}
// -----------------------------------------------------------------------------
double Utility::GetCVWeight(int type, double weightSplineTimesTune, double ppfx_cv, double nu_e, int nu_pdg, bool infv){

    // Always give weights of 1 to the data
    if (type == k_data || type == k_ext) return 1.0;

    double weight = 1.0;

    // Get the tune weight
    if (weight_tune) weight = weightSplineTimesTune;
    
    // Catch infinate/nan/unreasonably large tune weights
    CheckWeight(weight);

    // Get the PPFX CV flux correction weight
    double weight_flux = 1.0;
    if (weight_ppfx) weight_flux = ppfx_cv;

    CheckWeight(weight_flux);

    if (weight_ppfx) weight = weight * weight_flux;

    // Weight the below threshold events to zero. Current threhsold is 125 MeV
    if (type == k_mc && (nu_pdg == -12 || nu_pdg == 12) && nu_e <= 0.125) weight = 0.0;

    // This is the intrinsic nue weight that scales it to the standard overlay sample
    if (std::string(intrinsic_mode) == "intrinsic" && type == k_mc){
        
        // Kill off the out of fv events in the intrinsic nue sample (to avoid double counting)
        if (!infv) weight = 0.0; 
        else weight = weight * intrinsic_weight;
        
    }


    return weight;

}
// -----------------------------------------------------------------------------
void Utility::GetPiZeroWeight(double &weight, int pizero_mode, int nu_pdg, int ccnc, int npi0, double pi0_e){

    // Dont weight the nuecc events
    if ( (nu_pdg == 12 || nu_pdg == -12) && ccnc == k_CC) return;


    // Fix the normalisation
    if (pizero_mode == 1){
        
        if (npi0 > 0) {
            weight = weight * 0.759;
        }

    }
    // Try energy dependent scaling for pi0
    else if (pizero_mode == 2){
        
        if (npi0 > 0) {
            double pi0emax = 0.6;
            if (pi0_e > 0.1 && pi0_e < pi0emax){
                weight = weight * (1 - 0.4 * pi0_e);
            }
            else if (pi0_e > 0.1 && pi0_e >= pi0emax){
                weight = weight * (1 - 0.4 * pi0emax);
            }
            
        }
    }
    else {
        // Dont touch the weight
    }
    

}
// -----------------------------------------------------------------------------

bool Utility::GetDirectory(TFile* f, TDirectory* &d, TString string){
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
void Utility::CreateDirectory(std::string folder){

    std::string a = "if [ ! -d \"plots/";
    std::string b = "run" + std::string(run_period) + "/" + folder;
    std::string c = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
    std::string d = "run" + std::string(run_period) + "/" + folder;
    std::string e = "; fi";
    std::string command = a + b + c + d + e;
    system(command.c_str());
}
// -----------------------------------------------------------------------------
void Utility::Tabulate(bool inFV, std::string interaction, std::string classification, std::string pi0_classification, int type, std::vector<double> &counter_v, double weight) {

    if (type == k_mc){
        
        // Require in FV condition
        if (inFV){
            if (interaction == "nue_cc_qe")  counter_v.at(k_count_nue_cc_qe)  += weight; 
            if (interaction == "nue_cc_res") counter_v.at(k_count_nue_cc_res) += weight;
            if (interaction == "nue_cc_dis") counter_v.at(k_count_nue_cc_dis) += weight;
            if (interaction == "nue_cc_coh") counter_v.at(k_count_nue_cc_coh) += weight;
            if (interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_mec) += weight;

            if (interaction == "nue_bar_cc_qe")  counter_v.at(k_count_nuebar_cc_qe)  += weight; 
            if (interaction == "nue_bar_cc_res") counter_v.at(k_count_nuebar_cc_res) += weight;
            if (interaction == "nue_bar_cc_dis") counter_v.at(k_count_nuebar_cc_dis) += weight;
            if (interaction == "nue_bar_cc_coh") counter_v.at(k_count_nuebar_cc_coh) += weight;
            if (interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_mec) += weight;

            if (interaction == "numu_cc_qe"  )  counter_v.at(k_count_numu_cc_qe)  += weight;
            if (interaction == "numu_cc_res" )  counter_v.at(k_count_numu_cc_res) += weight;
            if (interaction == "numu_cc_dis" )  counter_v.at(k_count_numu_cc_dis) += weight;
            if (interaction == "numu_cc_coh" )  counter_v.at(k_count_numu_cc_coh) += weight;
            if (interaction == "numu_cc_mec" )  counter_v.at(k_count_numu_cc_mec) += weight;

            if (interaction == "numu_bar_cc_qe")  counter_v.at(k_count_numubar_cc_qe)  += weight;
            if (interaction == "numu_bar_cc_res") counter_v.at(k_count_numubar_cc_res) += weight;
            if (interaction == "numu_bar_cc_dis") counter_v.at(k_count_numubar_cc_dis) += weight;
            if (interaction == "numu_bar_cc_coh") counter_v.at(k_count_numubar_cc_coh) += weight;
            if (interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_mec) += weight;

            // These are all the nus, but now in the fv
            if (interaction == "nue_cc_qe" || interaction == "nue_cc_res" || interaction == "nue_cc_coh" || interaction == "nue_cc_dis" || interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_infv) += weight;
            if (interaction == "nue_bar_cc_qe" || interaction == "nue_bar_cc_res" || interaction == "nue_bar_cc_coh" || interaction == "nue_bar_cc_dis" || interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_infv) += weight;
            
            if (interaction == "numu_cc_qe" || interaction == "numu_cc_res" || interaction == "numu_cc_coh" || interaction == "numu_cc_dis" || interaction == "numu_cc_mec") counter_v.at(k_count_numu_cc_infv) += weight;
            if (interaction == "numu_bar_cc_qe" || interaction == "numu_bar_cc_res" || interaction == "numu_bar_cc_coh" || interaction == "numu_bar_cc_dis" || interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_infv) += weight;
        
            if (pi0_classification == "nue_cc")         counter_v.at(k_count_pi0_nue_cc_nopi0)    += weight;
            if (pi0_classification == "nue_cc_pi0")     counter_v.at(k_count_pi0_nue_cc_pi0)      += weight;
            if (pi0_classification == "nue_bar_cc")     counter_v.at(k_count_pi0_nuebar_cc_nopi0) += weight;
            if (pi0_classification == "nue_bar_cc_pi0") counter_v.at(k_count_pi0_nuebar_cc_pi0)   += weight;
            if (pi0_classification == "numu_cc")        counter_v.at(k_count_pi0_numu_cc_nopi0)   += weight;
            if (pi0_classification == "numu_cc_pi0")    counter_v.at(k_count_pi0_numu_cc_pi0)     += weight;
            if (pi0_classification == "nue_nc" || pi0_classification == "nue_bar_nc" || pi0_classification == "numu_nc" || pi0_classification == "numu_bar_nc")                 counter_v.at(k_count_pi0_nc_nopi0) += weight;
            if (pi0_classification == "nue_nc_pi0" || pi0_classification == "nue_bar_nc_pi0" || pi0_classification == "numu_nc_pi0" || pi0_classification == "numu_bar_nc_pi0") counter_v.at(k_count_pi0_nc_pi0) += weight;
        
        }

        // These are all the nus, but now in the cryostat volume
        if (interaction == "nue_cc_qe" || interaction == "nue_cc_res" || interaction == "nue_cc_coh" || interaction == "nue_cc_dis" || interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_incryo) += weight;
        if (interaction == "nue_bar_cc_qe" || interaction == "nue_bar_cc_res" || interaction == "nue_bar_cc_coh" || interaction == "nue_bar_cc_dis" || interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_incryo) += weight;
        
        if (interaction == "numu_cc_qe" || interaction == "numu_cc_res" || interaction == "numu_cc_coh" || interaction == "numu_cc_dis" || interaction == "numu_cc_mec") counter_v.at(k_count_numu_cc_incryo) += weight;
        if (interaction == "numu_bar_cc_qe" || interaction == "numu_bar_cc_res" || interaction == "numu_bar_cc_coh" || interaction == "numu_bar_cc_dis" || interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_incryo) += weight;
        
        // Classification
        if (classification == "nue_cc")           counter_v.at(k_count_nue_cc)           += weight;
        if (classification == "nuebar_cc")        counter_v.at(k_count_nuebar_cc)        += weight;
        if (classification == "nu_out_fv")        counter_v.at(k_count_nu_out_fv)        += weight;
        if (classification == "nc")               counter_v.at(k_count_nc)               += weight;
        if (classification == "nc_pi0")           counter_v.at(k_count_nc_pi0)           += weight;
        if (classification == "numu_cc")          counter_v.at(k_count_numu_cc)          += weight;
        if (classification == "numu_cc_pi0")      counter_v.at(k_count_numu_cc_pi0)      += weight;
        if (classification == "cosmic")           counter_v.at(k_count_cosmic)           += weight;
        if (classification == "unmatched")        counter_v.at(k_count_unmatched)        += weight;
        if (classification == "unmatched_nue")    counter_v.at(k_count_unmatched_nue)    += weight;
        if (classification == "cosmic_nue")       counter_v.at(k_count_cosmic_nue)       += weight;
        if (classification == "unmatched_nuebar") counter_v.at(k_count_unmatched_nuebar) += weight;
        if (classification == "cosmic_nuebar")    counter_v.at(k_count_cosmic_nuebar)    += weight;

        // Total selected MC events
        counter_v.at(k_count_total_mc) += weight;
    
    }
    else if (type == k_data) {
        counter_v.at(k_count_data) += weight; // The weight **should** always be 1 for these
    }
    else if (type == k_ext){
        counter_v.at(k_count_ext)  += weight;
    }
    else if (type == k_dirt){
        counter_v.at(k_count_dirt) += weight;
    }
    else {

        std::cout << "unkown type specified!!!  " << __PRETTY_FUNCTION__ << std::endl;
    }
    
}
// -----------------------------------------------------------------------------
double Utility::GetNuMIAngle(double px, double py, double pz, std::string direction){

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

    TVector3 beamdir = {0 , 0 , 1};;
    
    // Get the angle wrt to the beam
    if (direction == "beam") beamdir = {0 , 0 , 1};
    
    // Get the angle wrt to the target to detector direction
    else if (direction == "target") {
        beamdir = {5502, 7259, 67270};
        beamdir = beamdir.Unit(); // Get the direction
    }
    // NuMI acos(Pz/P)
    else if (direction == "numi_theta"){
        double p = std::sqrt(BeamCoords.X()*BeamCoords.X() +BeamCoords.Y()*BeamCoords.Y() +  BeamCoords.Z()*BeamCoords.Z() );
        return acos(BeamCoords.Z()/p) * (180 / 3.1415);
    }
    else if (direction == "numi_phi"){
        return atan2(BeamCoords.Y(), BeamCoords.X()) * 180 / 3.1415;
    } 
    else {
        std::cout << "Warning unknown angle type specified, you should check this" << std::endl;
    }
    
    double theta = BeamCoords.Angle(beamdir) * 180 / 3.1415926;


    // Create vectors to get the angle in the yz and xz planes
    TVector3 BeamCoords_yz = { 0, BeamCoords.Y(), BeamCoords.Z() }; // Angle upwards
    TVector3 BeamCoords_xz = { BeamCoords.X(), 0, BeamCoords.Z() }; // Angle across

    // if (theta > 50 ) std::cout <<"Theta: " << theta << "   UP: " << BeamCoords_yz.Angle(beam_dir) * 180 / 3.1415926 << "  Across: " << BeamCoords_xz.Angle(beam_dir) * 180 / 3.1415926 << std::endl;

    // std::cout << theta << std::endl;

    return theta;
}
// -----------------------------------------------------------------------------
bool Utility::in_fv(double x, double y, double z){
    
    // Require the vertex is in the boundary
    if ( x   >= config_v.at(k_config_x1) && x <= config_v.at(k_config_x2)
      && y   >= config_v.at(k_config_y1) && y <= config_v.at(k_config_y2)
      && z   >= config_v.at(k_config_z1) && z <= config_v.at(k_config_z2)  ){
        return true;   // pass
    }  
    else return false; // fail
}
// -----------------------------------------------------------------------------
void Utility::IncreaseLabelSize(TH1D *h, TCanvas *c){

    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);
}
// -----------------------------------------------------------------------------
void Utility::IncreaseLabelSize(TH2D *h, TCanvas *c){

    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.4);
    h->GetZaxis()->SetLabelSize(0.05);
    h->GetZaxis()->SetTitleSize(0.05);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.2);
    c->SetBottomMargin(0.13);
    h->SetMarkerSize(1.8);
    // gPad->SetGridx();
}
// -----------------------------------------------------------------------------
void Utility::Draw_Run_Period(TCanvas *c, double x1, double y1, double x2, double y2) {
    c->cd();

    //0.86, 0.915, 0.86, 0.915

    TPaveText *pt;

    if (std::string(run_period) == "1") {
        pt = new TPaveText(x1, y1, x2, y2, "NDC");
        pt->AddText("Run1");
        pt->SetTextColor(kRed + 2);
        pt->SetTextSize(0.04);
    }
    else if (std::string(run_period) == "3") {
        pt = new TPaveText(x1, y1, x2, y2, "NDC");
        pt->AddText("Run3");
        pt->SetTextColor(kBlue + 2);
    }
    else {
        pt = new TPaveText(x1, y1, x2, y2, "NDC");
        pt->AddText("RunXXX");
        pt->SetTextColor(kGreen + 2);
    }

    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.04);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void Utility::Draw_Data_MC_Ratio(TCanvas *c, double ratio, double x1, double y1, double x2, double y2){
    c->cd();

    // 0.34, 0.936, 0.34, 0.936

    TPaveText *pt;

    pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->AddText(Form("Data/MC Ratio: %2.2f", ratio));
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void Utility::Draw_Data_POT(TCanvas *c, double pot, double x1, double y1, double x2, double y2){
    c->cd();

    // 0.45, 0.915, 0.45, 0.915

    TPaveText *pt;

    // Change scale of POT
    double POT = pot / 1.0e20;

    pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->AddText(Form("MicroBooNE NuMI Data: %2.1f#times10^{20} POT", POT));
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void Utility::Draw_ubooneSim(TCanvas *c, double x1, double y1, double x2, double y2){
    c->cd();

    // 0.37, 0.92, 0.37, 0.92,

    TPaveText *pt;

    pt = new TPaveText(x1, y1, x2, y2,"NDC");
    pt->AddText("MicroBooNE Simulation");
    pt->SetTextColor(kBlack);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();

}
// -----------------------------------------------------------------------------
void Utility::SetTextProperties(TLatex* text){
    text->SetTextColor(kGray+2);
    text->SetNDC();
    text->SetTextSize(0.038);
    text->SetTextAlign(32);
}
// -----------------------------------------------------------------------------
void Utility::SetTPadOptions(TPad *topPad, TPad *bottomPad){
    topPad->SetBottomMargin(0.05);
    topPad->SetTopMargin(0.15);
    bottomPad->SetTopMargin(0.04);
    bottomPad->SetBottomMargin(0.25);
    bottomPad->SetGridy();
    topPad->SetLeftMargin(0.15);
    topPad->SetRightMargin(0.1);
    bottomPad->SetLeftMargin(0.15);
    bottomPad->SetRightMargin(0.1);
    topPad->Draw();
    bottomPad->Draw();
    topPad->cd();
}
// -----------------------------------------------------------------------------
void Utility::CheckPOT(){

    std::cout << "File: /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_overlay.root" << std::endl;

        // First we need to open the root file
        TFile * f = new TFile("/uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_overlay.root","READ");
        if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return; }

        // Get the tree
        TTree * mytree = (TTree*)f->Get("nuselection/SubRun");

        if (mytree == NULL){
            std::cout << "help can't get the branch so exiting..." << std::endl;
            gSystem->Exit(1);
        }

        float pot_sum = 0;
        float pot;
        mytree->SetBranchAddress("pot", &pot);

        // Loop over tree and get POT
        for (int i = 0; i < mytree->GetEntries(); i++) {
            mytree->GetEntry(i);
            //if (debug) std::cout << pot << std::endl;
            pot_sum = pot_sum + pot;
        }

        std::cout << "Total POT: " << pot_sum  << std::endl;

    // ----- check if the POT number in config.txt matched with pot_sum

    std::cout << std::endl;
    std::cout << "Checking the POT value: " << std::endl;

    std::ifstream config_file("../Analysis/config.txt"); // file with the POT values
    bool same_pot = false;

    std::string line; // saves the each line from config.txt
    
    while(getline(config_file,line)){
        
        if(line[0]=='R'){ // skips comment lines

            // splitting the line and taking the first arg as a TString and the second one as a float
            TString line_split(line);
            TObjArray *substrings = line_split.Tokenize(" ");
            float val = ((TObjString*)substrings->At(1))->GetString().Atof();
            TString label = ((TObjString*)substrings->At(0))->GetString();

            // calculate the ratio between new and old value, if =1 they are the same
            double ratio_pot = pot_sum/val;
            std::string ratio_pot_str = std::to_string(ratio_pot);

            // checks the initial label and if the ratio is one
            if(label == "Run1_MC_POT_CV" && ratio_pot_str == "1.000000") { // going back to string gives me a precision of 6 decimals
                std::cout << "The POT value of your file match the one in config.txt, continue running the code!" << std::endl;
                std::cout << std::endl;
            }

            // if POT values don't match, stop the code
            else if(label == "Run1_MC_POT_CV" && ratio_pot_str != "1.000000") {
                std::cout << std::endl;
                std::cout << "#################################################################" << std::endl;
                std::cout << " o    o     oo     ooo    oo   o   o   oo   o    ooooo " << std::endl;
                std::cout << " o    o    o  o    o  o   o o  o   o   o o  o   o      " << std::endl;
                std::cout << " o oo o   oooooo   o o    o  o o   o   o  o o   o    o " << std::endl;
                std::cout << " oo  oo   o    o   o  o   o   oo   o   o   oo    ooooo " << std::endl;
                std::cout << std::endl;
                std::cout << "POT value from your input file: " << pot_sum << std::endl;
                std::cout << "POT value saved in config.txt:  " << val << std::endl;
                std::cout << "The POT value of your file does not match the one in config.txt." << std::endl;
                std::cout << "Are you sure you want to continue?" << std::endl;
                std::cout << "#################################################################" << std::endl;
                std::cout << std::endl;
                //throw std::invalid_argument("The POT values don't match, please update config.txt");
            }
        }
    }

}
// -----------------------------------------------------------------------------
bool Utility::CheckHistogram(std::vector<std::string> vector, TString hist_name){
    
    // Add a temporary fix for inconsistent names
    if (std::string(hist_name) == "h_reco_vtx_x_sce" || std::string(hist_name) == "h_reco_vtx_y_sce" || std::string(hist_name) == "h_reco_vtx_z_sce")
        return true;


    for(unsigned int i=0; i<vector.size(); i++){

        if(vector.at(i).c_str() == hist_name) return true;

    }

    return false;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------