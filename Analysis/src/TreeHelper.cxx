#include "../include/TreeHelper.h"

// -----------------------------------------------------------------------------
void TreeHelper::Initialise(const char* run_period, const char * file_out, int weight_cfg ){

    std::cout << "Initalising Tree Helper..." << std::endl;

    std::string file_out_str = file_out;

    std::string file_name;

    // If the file name is empty then we use the default file name
    if (file_out_str == "empty") file_name = Form("files/nuexsec_selected_tree_mc_run%s.root", run_period);
    else file_name = "files/" + file_out_str;
    
    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
        f_nuexsec = new TFile( file_name.c_str(), "UPDATE");
    }
    
    // Set the weight settings
    if (weight_cfg == 0){
        weight_tune = false;
        weight_ppfx = false;
    }
    else if (weight_cfg == 1){
        weight_tune = true;
        weight_ppfx = true;
    }
    else if (weight_cfg == 2){
        weight_tune = true;
        weight_ppfx = false;
    }
    else if (weight_cfg == 3){
        weight_tune = false;
        weight_ppfx = true;
    }
    else {
        std::cout << "Unknown weight setting specified, using defaults" << std::endl;
    }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------