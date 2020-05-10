#include "../include/main.h"

int main(int argc, char *argv[]){

    auto start = std::chrono::high_resolution_clock::now(); // Start time of script

    bool using_slim_version        = false;
    bool make_histos               = false;
    bool run_selection             = true;
    bool area_norm                 = false;
    bool calc_cross_sec            = false;
    bool overwritePOT              = false; 
    bool run_sys                   = false;
    bool print                     = false;
    bool _print_mc                 = false;
    bool _print_data               = false;
    bool _print_ext                = false;
    bool _print_dirt               = false;

    // inputs 
    char * mc_file_name          = (char *)"empty";
    char * ext_file_name         = (char *)"empty";
    char * data_file_name        = (char *)"empty";
    char * dirt_file_name        = (char *)"empty";
    char * mc_file_name_out      = (char *)"empty";
    char * ext_file_name_out     = (char *)"empty";
    char * data_file_name_out    = (char *)"empty";
    char * dirt_file_name_out    = (char *)"empty";
    char * variation             = (char *)"empty";
    char * variation_file_name   = (char *)"empty";
    char * mc_tree_file_name_out = (char *)"empty";
    char * hist_file_name        = (char *)"empty";
    char * tree_file_name        = (char *)"empty";
    char * run_period            = (char *)"empty";
    int num_events{-1};
    int verbose{1}; // level 0 doesn't print cut summary, level 1 prints cut summary [default is 1 if unset]
    int weight{1};  // level 0 is no weights applied, level 1 (default) is all weights applied, level 2 is Genie Tune only, level 3 is PPFX CV only
    
    // Class instances
    xsecSelection::selection  _selection_instance;
    utility _utility;
    histogram_plotter _hplot;
    CrossSectionHelper _xsec;
    PrintHelper _phelper;
    // SystematicsHelper _systematics;

    std::string usage = "\nFirst run the selection with the options: \n\n\033[0;31m./nuexsec --run <run period num> [options (see below)]\033[0m \n\n"
    "\033[0;34m[--mc <mc file>]\033[0m                                       \033[0;32mThe input overlay root file\033[0m\n\n"
    "\033[0;34m[--data <data file>]\033[0m                                   \033[0;32mThe input on beam root file. Will also create a run subrun list in the ./files/ directory of selected data events.\033[0m\n\n"
    "\033[0;34m[--ext <ext file>]\033[0m                                     \033[0;32mThe input off beam root file\033[0m\n\n"
    "\033[0;34m[--dirt <dirt file>]\033[0m                                   \033[0;32mThe input dirt overlay root file\033[0m\n\n"
    "\033[0;33m[--mc_out <mc file output name>]\033[0m                       \033[0;32mThe output overlay root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--data_out <data file output name>]\033[0m                   \033[0;32mThe output on beam root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--ext_out <ext file output name>]\033[0m                     \033[0;32mThe output off beam root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--dirt_out <dirt file output name>]\033[0m                   \033[0;32mThe output dirt overlay root file name (will put in the ./files/ folder)\033[0m\n\n"
    "\033[0;33m[--var <variation file> <variation name>]\033[0m              \033[0;32m(first arg) Path to variation file name, (second arg) the variation name -- overwrites the default MC option\033[0m\n\n"
    "\033[0;33m[-n <num events>]\033[0m                                      \033[0;32mThe number of events to run over. This is good for checking if the code doesn't segfault. All the POT scalings will not work.\033[0m\n\n"
    "\033[0;33m[--weight <weight setting>]\033[0m                            \033[0;32mChange the Weight level. level 0 is no weights applied, level 1 (default) is all weights applied, level 2 is Genie Tune only, level 3 is PPFX CV only \033[0m\n\n"
    "\033[0;33m[--slim]\033[0m                                               \033[0;32mWhen this extension is added, the histogram helper class is not initalised and no histograms will be filled or saved. This is to speed up the selection code if you just want to run the selection.\033[0m\n\n"
    "\033[0;33m[--verbose <verbose level>]\033[0m                            \033[0;32m0 does not print the selection cut results, 1 (default) currently prints everything\033[0m\n\n"
    "-------------------------------------------------------"
    "\n\nTo make the histograms after running selection, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --hist <input merged nuexsec file> [options (see below)]\033[0m \n\n"
    "\033[0;33m[--weight <weight setting>]\033[0m                            \033[0;32mChange the Weight level to dislay on the plots. Should be used in conjunction with the setting used in the selection stage. level 0 is no weights applied, level 1 (default) is all weights applied, level 2 is Genie Tune only, level 3 is PPFX CV only \033[0m\n\n"
    "\033[0;33m[--area]\033[0m                                               \033[0;32mArea normalise all the histograms\033[0m\n\n"
    "\033[0;33m[--var dummy <variation name>]\033[0m              \033[0;32m(first arg) this argument is already input from the hist option, use something like -dummy- as a placeholder, (second arg) the variation name\033[0m\n\n"
    "The <input merged nuexsec file> corresponds to hadd merged file of the mc, data, ext and dirt. See the bash script merge_run1_files.sh for more details\n\n"
    "-------------------------------------------------------"
    "\n\nTo run the cross section calculation code, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --xsec <input merged nuexsec tree file> [options (see below)]\033[0m \n\n"
    "The <input merged nuexsec ttree file> corresponds to merged ttree file of the mc, data, ext and dirt. See the bash script merge_uneaventrees.C for more details\n\n";

    std::string usage2 = "\n\nTo run the detector systematics code, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --sys\033[0m \n\n"
    "This will run the detector systematics plotting code\n\n";


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

        // Histogram Mode
        if (strcmp(arg, "--hist") == 0) {
            std::cout << "Making Histograms, file to make histograms with: "<< argv[i+1] << std::endl;
            make_histos = true;
            run_selection = false; // switch this bool out
            hist_file_name = argv[i+1];
        }

        // XSec mode
        if (strcmp(arg, "--xsec") == 0) {
            std::cout << "Calculating the Cross Section with input file name: " << argv[i+1] << std::endl;
            calc_cross_sec = true;
            run_selection = false; // switch this bool out
            tree_file_name = argv[i+1];
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

        // Weight Settings
        if (strcmp(arg, "--weight") == 0){
            std::cout << "Running with weight configuration setting of: " << argv[i+1] << std::endl;
            weight = atoi(argv[i+1]);
        }

        // Whats the verbose?
        if (strcmp(arg, "-v") == 0 || strcmp(arg, "--verbose") == 0 || strcmp(arg, "--v") == 0){
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

        // Help!
        if (strcmp(arg, "--h") == 0 || strcmp(arg, "-h") == 0|| strcmp(arg, "--help") == 0 || strcmp(arg, "--usage") == 0){
            std::cout << usage << std::endl; 
            exit(1);
        }

        // Area Normalise the histograms
        if (strcmp(arg, "--area") == 0) {
            std::cout << "Area Normalising the histograms"<< std::endl;
            area_norm = true;
        }

        // Variation file
        if (strcmp(arg, "--var") == 0){
            std::cout << "Using Variation: " << argv[i+2] << std::endl;
            
            variation = argv[i+2];
            mc_file_name = argv[i+1];
            
            mc_file_name_out      = Form("nuexsec_mc_run%s_%s.root", run_period, variation);
            mc_tree_file_name_out = Form("nuexsec_selected_tree_mc_run%s_%s.root", run_period, variation);
            
            std::cout  << "Output filename will be overwritten with name: "      << mc_file_name_out << std::endl;
            std::cout << "Output tree filename will be overwritten with name: " << mc_tree_file_name_out << std::endl;
            
            overwritePOT = true;
        }

        // Systematics
        if (strcmp(arg, "--sys") == 0){
            // std::cout << "Using Systematics plotting code" << std::endl;
            
            run_sys = true;
            run_selection = false;
        }

        // Only run the print function
        if (strcmp(arg, "--printonly") == 0) {
            run_selection = false;
            print = true;
        }

        // Print all
        if (strcmp(arg, "--printall") == 0) {
            print = true;
            _print_mc                 = true;
            _print_data               = true;
            _print_ext                = true;
            _print_dirt               = true;
        }
        
        // Print MC
        if (strcmp(arg, "--printmc") == 0) {
            print = true;
            _print_mc = true;
        }
        
        // Print Data
        if (strcmp(arg, "--printdata") == 0) {
            print = true;
            _print_data = true;
        }

        // Print ext
        if (strcmp(arg, "--printext") == 0) {
            print = true;
            _print_ext = true;
        }

        // Print Dirt
        if (strcmp(arg, "--printdirt") == 0) {
            print = true;
            _print_dirt = true;
        }
   
    }

    // Add catches for default input
    if ((run_period == "empty") ){
        std::cout << "\nError, must provide a run period as input!\n" << std::endl;
        std::cout << "USAGE:" << usage << std::endl; 
        exit(1);
    }

    // -------------------------------------------------------------------------
    
    // Configure the utility class
    _utility.Initalise(variation, overwritePOT, run_period);

    // -------------------------------------------------------------------------

    // Initialise the selction script
    if (run_selection) _selection_instance.Initialise(mc_file_name, ext_file_name, data_file_name, dirt_file_name,
                                                      mc_file_name_out, ext_file_name_out, data_file_name_out, dirt_file_name_out,
                                                      mc_tree_file_name_out, _utility, using_slim_version, num_events, run_period, verbose, weight );

    
    // Print the selection results
    if (print) _phelper.Initialise(run_period, "empty", weight, _print_mc, _print_data, _print_ext, _print_dirt, _utility );

    // Run the make histogram function
    if (make_histos) _hplot.MakeHistograms(hist_file_name, run_period, weight, area_norm, _utility, variation);

    // Run the calculate cross section function
    if (calc_cross_sec) _xsec.Initialise(run_period, tree_file_name, _utility);

    // -------------------------------------------------------------------------
    // Finished!
    std::cout << "\033[0;32m*** \t Exiting C++ Code... \t *** \033[0m" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();  // end time of script
    auto duration_sec = std::chrono::duration_cast<std::chrono::seconds>(stop - start); // time taken to run script
    auto duration_min = std::chrono::duration_cast<std::chrono::minutes>(stop - start); // time taken to run script
    std::cout << "Time taken by function: " << duration_sec.count() << " seconds" << std::endl; 
    std::cout << "Time taken by function: " << duration_min.count() << " minutes" << std::endl; 
    
    exit(0);
    return 0;
}
