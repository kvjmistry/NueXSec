#include "../include/main.h"

int main(int argc, char *argv[]){

    auto start = std::chrono::high_resolution_clock::now(); // Start time of script

    const char * input_config_file_name;
    bool using_default_config      = true;
    bool using_slim_version        = false;
    bool make_histos               = false;

    char * mc_file_name          = (char *)"empty";
    char * ext_file_name         = (char *)"empty";
    char * data_file_name        = (char *)"empty";
    char * dirt_file_name        = (char *)"empty";
    char * variation_file_name   = (char *)"empty";

    // -------------------------------------------------------------------------
    // Loop over input arguments
    for (int i =1; i < argc; i++){
        auto const arg = argv[i];
        //std::cout << arg << std::endl; //this is for debugging
        
        // Slim input
        if(strcmp(arg, "--slim") == 0) {
            std::cout << "Running with slim mode"<< std::endl;
            using_slim_version = true;
            std::cout << " *** \t Running with Slimmed Selection (no histograms will be made)\t *** " << std::endl;
        }

        if(strcmp(arg, "--hist") == 0) {
            std::cout << "Making Histograms"<< std::endl;
            make_histos = true;
        }
        
        // Configuration file
        if(strcmp(arg, "-c") == 0){
            std::cout << "Running with input config file: " << argv[i+1]<< std::endl;
            using_default_config = false;
            input_config_file_name = argv[i+1];
        }
        
        // MC file
        if(strcmp(arg, "--mc") == 0) {
            std::cout << "Running with MC file: " << argv[i+1]<< std::endl;
            mc_file_name = argv[i+1];
        }

        // EXT file
        if(strcmp(arg, "--ext") == 0){
            std::cout << "Running with EXT file: " << argv[i+1]<< std::endl;
            ext_file_name = argv[i+1];
        }

        // Data file
        if(strcmp(arg, "--data") == 0){
            std::cout << "Running with Data file: " << argv[i+1]<< std::endl;
            data_file_name = argv[i+1];
        }

        // Dirt file overlay or not?
        if(strcmp(arg, "--dirt") == 0){
            std::cout << "Running with Dirt file: " << argv[i+1]<< std::endl;
            dirt_file_name = argv[i+1];
        }

        // Variation file
        if(strcmp(arg, "--var") == 0){
            std::cout << "Running with Systematic Variation file: " << argv[i+1]<< std::endl;
            variation_file_name = argv[i+1];
            std::string variation_type = variation_file_name;
        }
        if (strcmp(arg, "--h") == 0 || strcmp(arg, "-h") == 0|| strcmp(arg, "--help") == 0 || strcmp(arg, "--usage") == 0){
            std::cout << " \n ./nuexsec --mc <mc file> [--data <data file>] [--ext <ext file>] [--dirt <dirt file>] [--var <variation file>] [-c <input config file>] [--slim] [--hist] \n " << std::endl; 
            exit(1);
        }
        
    }
    // if(argc < 2 )  { std::cout << " \n Please include the input file path \n " << std::endl; exit(1); }

    // -------------------------------------------------------------------------
    std::vector<double> config;
    std::vector<double> input_config;
    xsecSelection::selection  _selection_instance;

    // Configure the cut values
    utilityNS::utility _utility;
    std::vector<double> default_config = _utility.configure_cuts(
            _x1, _x2, _y1, _y2, _z1, _z2,
            flash_pe_threshold,
            flash_time_start,
            flash_time_end,
            tolerance,
            shwr_nue_tolerance,
            trk_nue_tolerance,
            shwr_hit_threshold,
            shwr_hit_threshold_collection,
            tolerance_open_angle_min,
            tolerance_open_angle_max,
            tolerance_dedx_min,
            tolerance_dedx_max,
            dist_tolerance,
            pfp_hits_length_tolerance,
            ratio_tolerance,
            detector_variations
            );
    // -------------------------------------------------------------------------

    // Check if a config file has been used for the cuts
    std::ifstream input_config_file;
    if(using_default_config == false) {
        input_config_file.open(input_config_file_name);
        
        // Bad file
        if( !input_config_file.is_open()) {
            std::cout << "*** \t File did not open! Quitting \t ***" << std::endl;
            exit(1);
        }
        
        // Good file
        if (input_config_file.is_open()) {
            std::cout << "*** \t File Opened \t *** " << std::endl;
            std::string line;

            // Loop over the lines in the file
            while (!input_config_file.eof()) {
                std::getline (input_config_file, line);
                if(!line.empty()) {
                    input_config.push_back(std::stof(line));
                }
            }
            input_config_file.close();
            std::cout << "*** \t File Closed \t *** " << std::endl;
        }
    }
    // -------------------------------------------------------------------------
    // Set the config
    if(using_default_config == true) {
        config = default_config;
        std::cout << "(Default) Parameter Config from header" << std::endl;
    }
    
    if(using_default_config == false) {
        config = input_config;
        std::cout << "Paramerter Config from input file" << std::endl;
    }
    // -------------------------------------------------------------------------

    // Initialise the selction script
    _selection_instance.xsecSelection::selection::Initialise(mc_file_name, ext_file_name, data_file_name, dirt_file_name, variation_file_name, config, using_slim_version );

    // now save all the outputs to file
    if (!using_slim_version) _selection_instance.xsecSelection::selection::SavetoFile();
    
    if (make_histos) _selection_instance.xsecSelection::selection::MakeHistograms();

    // -------------------------------------------------------------------------
    // Finished!
    std::cout << "\033[0;32m*** \t Exiting C++ Code... \t *** \033[0m" << std::endl;
    std::cout << "\033[0;31m*** \t Warning shift of 1us is still set in flash times \t *** \033[0m" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();  // end time of script
	auto duration_sec = std::chrono::duration_cast<std::chrono::seconds>(stop - start); // time taken to run script
	auto duration_min = std::chrono::duration_cast<std::chrono::minutes>(stop - start); // time taken to run script
	std::cout << "Time taken by function: " << duration_sec.count() << " seconds" << std::endl; 
	std::cout << "Time taken by function: " << duration_min.count() << " minutes" << std::endl; 
    
    return 0;
}
