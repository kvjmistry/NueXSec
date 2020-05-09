#include "../include/SystematicsHelper.h"

// -----------------------------------------------------------------------------
void SystematicsHelper::Initialise(const char *_run_period, utility _utility){

    std::cout << "Initalising Systematics Helper..." << std::endl;
    _util = _utility;

    // Get the POT of the variations from the file
    GetPOT(_run_period);

    // Set the run period
    run_period = std::string(_run_period);

    // Get the variation files
    for (unsigned int l =0; l < var_string.size(); l++){
        f_vars.at(l) = new TFile( Form("files/nuexsec_run%s_%s_merged", _run_period, var_string.c_str() ), "READ");

        // Maybe add something here to pickup non processed variation
    }

    // Now loop over events and caluclate the cross section
    MakeHistograms();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::MakeHistograms(){






}
// -----------------------------------------------------------------------------
void SystematicsHelper::GetPOT(const char* run_period){

    std::cout << "Getting the POT for the variation files" << std::endl;

    POT_v.resize(k_vars_MAX, 1.0);

    std::string line;

    std::string varname;
    std::string value;

    std::string POT_run_config = "Run" + std::string(run_period) + "MC_POT_";
    
    // Loop over the config ist
    for (unsigned int p = 0; p < var_string.size(); p++){

        std::ifstream myfile ("config.txt");

        if (myfile.is_open()) {
            
            // std::cout << var_string.at(p) <<  std::endl;

            // Loop over lines in file
            while ( getline (myfile,line) ) {

                std::istringstream ss(line);
                ss >> varname >> value;

                // Found the correct variation file 
                std::string POT_name_match = POT_run_config + var_string.at(p);
                if (varname == POT_name_match ) {
                    std::cout << "Found match for: " << varname << " "<< std::stod(value) <<  std::endl;
                    POT_v.at(p)= std::stod(value);
                    break;
                }
                
            }
           
        }
        else std::cout << "Unable to open config file, bad things are going to happen..." << std::endl; 

        myfile.close();
    }

}
// -----------------------------------------------------------------------------