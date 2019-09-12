#include "../include/selection.h"

namespace xsecSelection {

// Main function for selection
void selection::Initialise( const char * mc_file,
                            const char * ext_file,
                            const char * data_file,
                            const char * dirt_file,
                            const char * variation_file,
                            const std::vector<double> config
                            ){
    
    // Print the input files
    std::cout <<
    "MC   File Path:      " << mc_file        <<"\n" <<
    "Ext  File Path:      " << ext_file       <<"\n" <<
    "Data File Path:      " << data_file      <<"\n" <<
    "Dirt File Path:      " << dirt_file      <<"\n" <<
    "Variation File Path: " << variation_file <<"\n" <<
    std::endl;

    // Now get the files, if file isnt specified then set bool to skip
    bool_use_mc        = _utility_instance.GetFile(f_mc,        mc_file);
    bool_use_ext       = _utility_instance.GetFile(f_ext,       ext_file);
    bool_use_data      = _utility_instance.GetFile(f_data,      data_file);
    bool_use_dirt      = _utility_instance.GetFile(f_dirt,      dirt_file);
    bool_use_variation = _utility_instance.GetFile(f_variation, variation_file);


    // Get MC variables
    if (bool_use_mc){
        // utility::GetTree(f_mc, mytree, "AnalyzeTPCO/tree");
        // utility::GetTree(f_mc, optree, "AnalyzeTPCO/optical_tree");
        // utility::GetTree(f_mc, mctree, "AnalyzeTPCO/mcparticletree");
        // utility::GetTree(f_mc, mctruth_counter_tree, "AnalyzeTPCO/mctruth_counter_tree");

    } // End getting MC variables
    

} // END make_selection function
} // END NAMESPACE xsecSelection

