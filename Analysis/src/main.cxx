#include "../include/main.h"

int main(int argc, char *argv[]){

    auto start = std::chrono::high_resolution_clock::now(); // Start time of script

    bool using_slim_version        = false;
    bool make_histos               = false;

    // inputs 
    char * mc_file_name          = (char *)"empty";
    char * ext_file_name         = (char *)"empty";
    char * data_file_name        = (char *)"empty";
    char * dirt_file_name        = (char *)"empty";
    char * variation_file_name   = (char *)"empty";
    char * hist_file_name        = (char *)"empty";
    char * run_period            = (char *)"empty";
    int num_events{-1};
    int verbose{1}; // level 0 doesn't print cut summary, level 1 prints cut summary [default is 1 if unset]

    // Configurations from main.h will be stored in here and passed to the selection
    std::vector<double> config;
    
    // Class instances
    xsecSelection::selection  _selection_instance;
    utility _utility;

    std::string usage = "\n First run the selection with the options: \n\n ./nuexsec --run <run period num> [--mc <mc file>] [--data <data file>] [--ext <ext file>] [--dirt <dirt file>] [--var <variation file>] [-n <num events>] [--slim] [--verbose <verbose level>] \n\n After this, to make the histograms after running selection, run: \n\n ./nuexsec --run <run period num> --hist <input merged nuexsec file> \n ";

    // -------------------------------------------------------------------------
    // Loop over input arguments
    for (int i =1; i < argc; i++){
        auto const arg = argv[i];
        //std::cout << arg << std::endl; // This is for debugging
        
        // Slim input
        if (strcmp(arg, "--slim") == 0) {
            std::cout << "Running with slim mode"<< std::endl;
            using_slim_version = true;
            std::cout << " *** \t Running with Slimmed Selection (no histograms will be made)\t *** " << std::endl;
        }

        if (strcmp(arg, "--hist") == 0) {
            std::cout << "Making Histograms, file to make histograms with: "<< argv[i+1] << std::endl;
            make_histos = true;
            hist_file_name = argv[i+1];
        }
        
        // MC file
        if (strcmp(arg, "--mc") == 0) {
            std::cout << "Running with MC file: " << argv[i+1] << std::endl;
            mc_file_name = argv[i+1];
        }

        // EXT file
        if (strcmp(arg, "--ext") == 0){
            std::cout << "Running with EXT file: " << argv[i+1] << std::endl;
            ext_file_name = argv[i+1];
        }

        // Data file
        if (strcmp(arg, "--data") == 0){
            std::cout << "Running with Data file: " << argv[i+1] << std::endl;
            data_file_name = argv[i+1];
        }

        // Dirt file overlay or not?
        if (strcmp(arg, "--dirt") == 0){
            std::cout << "Running with Dirt file: " << argv[i+1] << std::endl;
            dirt_file_name = argv[i+1];
        }

        // Variation file
        if (strcmp(arg, "--var") == 0){
            std::cout << "Running with Systematic Variation file: " << argv[i+1] << std::endl;
            variation_file_name = argv[i+1];
            std::string variation_type = variation_file_name;
        }

        // Whats the verbose?
        if (strcmp(arg, "-v") == 0 || strcmp(arg, "--verbose") == 0){
            std::cout << "Setting Verbose Level to : " << argv[i+1] << std::endl;
            verbose = atoi(argv[i+1]);
        }

        // Max number of events specified?
        if (strcmp(arg, "-n") == 0 || strcmp(arg, "--n") == 0){
            std::cout << "Running with a maximum of : " << argv[i+1] << " events" <<std::endl;
            num_events = atoi(argv[i+1]);
        }

        // Set the run period
        if (strcmp(arg, "--run") == 0){
            std::cout << "Setting the run period as : run" << argv[i+1] <<std::endl;
            run_period = argv[i+1];
        }

        if (strcmp(arg, "--h") == 0 || strcmp(arg, "-h") == 0|| strcmp(arg, "--help") == 0 || strcmp(arg, "--usage") == 0){
            std::cout << usage << std::endl; 
            exit(1);
        }
        
    }

    // Add catches for default input
    if ((run_period == "empty") ){
        std::cout << "Error, must provide a run period as input!" << std::endl;
        std::cout << "USAGE:" << usage << std::endl; 
        exit(1);
    }

    // -------------------------------------------------------------------------
    
    // Configure the cut values
    config = _utility.configure(
        _Run1_MC_POT,
        _Run1_Dirt_POT,
        _Run1_Data_POT,
        _Run1_Data_trig,
        _Run1_EXT_trig,
        _x1, _x2, _y1, _y2, _z1, _z2
        );

    // -------------------------------------------------------------------------

    // Initialise the selction script
    if (!make_histos) _selection_instance.xsecSelection::selection::Initialise(mc_file_name, ext_file_name, data_file_name, dirt_file_name, variation_file_name, config, using_slim_version, num_events, run_period, verbose );
    else _selection_instance.xsecSelection::selection::MakeHistograms(hist_file_name, run_period, config);

    // -------------------------------------------------------------------------
    // Finished!
    std::cout << "\033[0;32m*** \t Exiting C++ Code... \t *** \033[0m" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();  // end time of script
	auto duration_sec = std::chrono::duration_cast<std::chrono::seconds>(stop - start); // time taken to run script
	auto duration_min = std::chrono::duration_cast<std::chrono::minutes>(stop - start); // time taken to run script
	std::cout << "Time taken by function: " << duration_sec.count() << " seconds" << std::endl; 
	std::cout << "Time taken by function: " << duration_min.count() << " minutes" << std::endl; 
    
    return 0;
}
