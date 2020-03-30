#include "../include/main.h"

int main(int argc, char *argv[]){

    auto start = std::chrono::high_resolution_clock::now(); // Start time of script

    bool using_slim_version        = false;
    bool make_histos               = false;
    bool run_selection             = true;

    // inputs 
    char * mc_file_name          = (char *)"empty";
    char * ext_file_name         = (char *)"empty";
    char * data_file_name        = (char *)"empty";
    char * dirt_file_name        = (char *)"empty";
    char * mc_file_name_out      = (char *)"empty";
    char * ext_file_name_out     = (char *)"empty";
    char * data_file_name_out    = (char *)"empty";
    char * dirt_file_name_out    = (char *)"empty";
    char * variation_file_name   = (char *)"empty";
    char * hist_file_name        = (char *)"empty";
    char * run_period            = (char *)"empty";
    int num_events{-1};
    int verbose{1}; // level 0 doesn't print cut summary, level 1 prints cut summary [default is 1 if unset]
    int weight{1};  // level 0 is no weights applied, level 1 (default) is all weights applied, level 2 is Genie Tune only, level 3 is PPFX CV only

    // Configurations from main.h will be stored in here and passed to the selection
    std::vector<double> config;
    
    // Class instances
    xsecSelection::selection  _selection_instance;
    utility _utility;
    histogram_plotter _hplot;

    std::string usage = "\nFirst run the selection with the options: \n\n\033[0;31m./nuexsec --run <run period num> [options (see below)]\033[0m \n\n"
    "\033[0;34m[--mc <mc file>]\033[0m                                       \033[0;32mThe input overlay root file\033[0m\n\n"
    "\033[0;34m[--data <data file>]\033[0m                                   \033[0;32mThe input on beam root file. Will also create a run subrun list in the ./files/ directory of selected data events.\033[0m\n\n"
    "\033[0;34m[--ext <ext file>]\033[0m                                     \033[0;32mThe input off beam root file\033[0m\n\n"
    "\033[0;34m[--dirt <dirt file>]\033[0m                                   \033[0;32mThe input dirt overlay root file\033[0m\n\n"
    "\033[0;33m[--mc_out <mc file output name>]\033[0m                       \033[0;32mThe output overlay root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--data_out <data file output name>]\033[0m                   \033[0;32mThe output on beam root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--ext_out <ext file output name>]\033[0m                     \033[0;32mThe output off beam root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--dirt_out <dirt file output name>]\033[0m                   \033[0;32mThe output dirt overlay root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--var <variation file>]\033[0m                               \033[0;31mNot yet supported\033[0m\n\n"
    "\033[0;33m[-n <num events>]\033[0m                                      \033[0;32mThe number of events to run over. This is good for checking if the code doesn't segfault. All the POT scalings will not work.\033[0m\n\n"
    "\033[0;33m[--weight <weight setting>]\033[0m                            \033[0;32mChange the Weight level. level 0 is no weights applied, level 1 (default) is all weights applied, level 2 is Genie Tune only, level 3 is PPFX CV only \033[0m\n\n"
    "\033[0;33m[--slim]\033[0m                                               \033[0;32mWhen this extension is added, the histogram helper class is not initalised and no histograms will be filled or saved. This is to speed up the selection code if you just want to run the selection.\033[0m\n\n"
    "\033[0;33m[--verbose <verbose level>]\033[0m                            \033[0;32m0 does not print the selection cut results, 1 (default) currently prints everything\033[0m\n\n"
    "\n\nTo make the histograms after running selection, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --hist <input merged nuexsec file>\033[0m \n\n"
    "The <input merged nuexsec file> corresponds to hadd merged file of the mc, data, ext and dirt. See the bash script merge_run1_files.sh for more details\n\n";

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
            run_selection = false; // switch this bool out
            hist_file_name = argv[i+1];
        }
        
        // MC file
        if (strcmp(arg, "--mc") == 0) {
            std::cout << "Running with MC file: " << argv[i+1] << std::endl;
            mc_file_name = argv[i+1];
        }

        // Overwrite output mc file name
        if (strcmp(arg, "--mc_out") == 0) {
            std::cout << "New Output MC File name: " << argv[i+1] << std::endl;
            mc_file_name_out = argv[i+1];
        }

        // EXT file
        if (strcmp(arg, "--ext") == 0){
            std::cout << "Running with EXT file: " << argv[i+1] << std::endl;
            ext_file_name = argv[i+1];
        }

        // Overwrite output ext file name
        if (strcmp(arg, "--ext_out") == 0) {
            std::cout << "New Output EXT File name: " << argv[i+1] << std::endl;
            ext_file_name_out = argv[i+1];
        }

        // Data file
        if (strcmp(arg, "--data") == 0){
            std::cout << "Running with Data file: " << argv[i+1] << std::endl;
            data_file_name = argv[i+1];
        }
        
        // Overwrite output data file name
        if (strcmp(arg, "--data_out") == 0) {
            std::cout << "New Output Data File name: " << argv[i+1] << std::endl;
            data_file_name_out = argv[i+1];
        }

        // Dirt file overlay
        if (strcmp(arg, "--dirt") == 0){
            std::cout << "Running with Dirt file: " << argv[i+1] << std::endl;
            dirt_file_name = argv[i+1];
        }

        // Overwrite output dirt file name
        if (strcmp(arg, "--dirt_out") == 0) {
            std::cout << "New Output Dirt File name: " << argv[i+1] << std::endl;
            dirt_file_name_out = argv[i+1];
        }

        // Variation file
        if (strcmp(arg, "--var") == 0){
            std::cout << "Running with Systematic Variation file: " << argv[i+1] << std::endl;
            variation_file_name = argv[i+1];
            std::string variation_type = variation_file_name;
        }

        // Variation file
        if (strcmp(arg, "--weight") == 0){
            std::cout << "Running with weight configuration setting of: " << argv[i+1] << std::endl;
            weight = atoi(argv[i+1]);
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
    if (run_selection) _selection_instance.xsecSelection::selection::Initialise(mc_file_name, ext_file_name, data_file_name, dirt_file_name, mc_file_name_out, ext_file_name_out, data_file_name_out, dirt_file_name_out, variation_file_name, config, using_slim_version, num_events, run_period, verbose, weight );
    
    // Run the make histogram function
    if (make_histos) _hplot.MakeHistograms(hist_file_name, run_period, config, weight);

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
