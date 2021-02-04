#include "../include/Main.h"

int main(int argc, char *argv[]){

    auto start = std::chrono::high_resolution_clock::now(); // Start time of script

    // Class instances
    xsecSelection::Selection  _selection_instance;
    Utility              _utility;
    HistogramPlotter     _hplot;
    CrossSectionHelper   _xsec;
    PrintHelper          _phelper;
    SystematicsHelper    _syshelper;
    UtilityPlotter       _uplot;

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
    "\033[0;33m[--weight_tune <weight mode>]\033[0m                          \033[0;32mTurn on/off the GENIE Tune 1 == on, 0 == off \033[0m\n\n"
    "\033[0;33m[--weight_ppfx <weight mode>]\033[0m                          \033[0;32mTurn on/off the PPFX CV 1 == on, 0 == off \033[0m\n\n"
    "\033[0;33m[--weight_dirt <weight mode>]\033[0m                          \033[0;32mTurn on/off the weighting of dirt 1 == on, 0 == off \033[0m\n\n"
    "\033[0;33m[--weight_ext <weight mode>]\033[0m                           \033[0;32mTurn on/off the weighting of ext 1 == on, 0 == off \033[0m\n\n"
    "\033[0;33m[--weight_pi0 <weight mode>]\033[0m                           \033[0;32mTurn on/off the weighting of pi0 0 == off, 1 == norm factor, 2 == E dep. scaling \033[0m\n\n"
    "\033[0;33m[--notuning]\033[0m                                           \033[0;32mTurn off ppfx, genie tune and pi0 tunings. \033[0m\n\n"
    "\033[0;33m[--slim]\033[0m                                               \033[0;32mWhen this extension is added, the histogram helper class is not initalised and no histograms will be filled or saved. This is to speed up the selection code if you just want to run the selection.\033[0m\n\n"
    "\033[0;33m[--intrinsic <intrinsic mode>]\033[0m                         \033[0;32mWhen this is on, we overwrite the signal events that get written to file with events from an intrinsic file. Events are also weighted by an additional intrinsic factor. Must be used with --mc option. Options are intrinsic.\033[0m\n\n"
    "\033[0;33m[--tunemec]\033[0m                                            \033[0;32mApply a 1.5 scale factor to MEC events. Default is off \033[0m\n\n"
    // "\033[0;33m[--verbose <verbose level>]\033[0m                            \033[0;32mDoes not print the selection cut results, 1 (default) currently prints everything\033[0m\n\n"
    "-------------------------------------------------------"
    "\n\nTo print the results of the selection, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> [options (see below)]\033[0m \n\n"
    "\033[0;33m[--printmc]\033[0m                                            \033[0;32mPrints the selection results if MC ran \033[0m\n\n"
    "\033[0;33m[--printdata]\033[0m                                          \033[0;32mPrints the selection results if Data ran \033[0m\n\n"
    "\033[0;33m[--printext]\033[0m                                           \033[0;32mPrints the selection results if EXT ran \033[0m\n\n"
    "\033[0;33m[--printdirt]\033[0m                                          \033[0;32mPrints the selection results if Dirt ran \033[0m\n\n"
    "\033[0;33m[--printall]\033[0m                                           \033[0;32mPrints the selection results if mc, data, ext, dirt ran \033[0m\n\n"
    "\033[0;34m[--mc <mc file>]\033[0m                                       \033[0;32m\033[0mEntirly optional, but allows you to override the MC file that gets printed\n\n"
    "\033[0;33m[--printonly]\033[0m                                          \033[0;32mTurns off running of the selection, turns on print function. Good for checking selection results and not having to run the selection again \033[0m\n\n";
    
    std::string usage2 = "-------------------------------------------------------"
    "\n\nTo make the histograms after running selection, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --hist <input merged nuexsec file> [options (see below)]\033[0m \n\n"
    // "\033[0;33m[--weight <weight setting>]\033[0m                            \033[0;32mChange the Weight level to dislay on the plots. Should be used in conjunction with the setting used in the selection stage. level 0 is no weights applied, level 1 (default) is all weights applied, level 2 is Genie Tune only, level 3 is PPFX CV only \033[0m\n\n"
    "\033[0;33m[--area]\033[0m                                               \033[0;32mArea normalise all the histograms\033[0m\n\n"
    "\033[0;33m[--var dummy <variation name>]\033[0m                         \033[0;32m(first arg) this argument is already input from the hist option, use something like -dummy- as a placeholder, (second arg) the variation name\033[0m\n\n"
    "The <input merged nuexsec file> corresponds to hadd merged file of the mc, data, ext and dirt. See the bash script merge_run1_files.sh for more details\n\n"
    "-------------------------------------------------------"
    "\n\nTo run the cross section calculation code, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --xsec <input merged nuexsec tree file> [options (see below)]\033[0m \n\n"
    "The <input merged nuexsec ttree file> corresponds to merged ttree file of the mc, data, ext and dirt. See the bash script merge_uneaventrees.C for more details\n\n"
    "\033[0;34m[--xsecmode <cross-section mode> ]\033[0m                     \033[0;32mThe input mode of xsec code to run. Options are default or reweight. Labels are optional, but can be all/unisim/ppfx/genie/reint \033[0m\n\n"
    "\033[0;34m[--xseclabel <cross-section label> ]\033[0m                   \033[0;32mThe systematic variation to run over in the xsec calc. Default is all. Options are all/unisim/ppfx/genie/reint \033[0m\n\n"
    "\033[0;34m[--xsecplot <cross-section plot mode> ]\033[0m                \033[0;32mChoose to use the cross section code to reweight all the cut plots [very slow]. Options are default or rw_cuts. \033[0m\n\n"
    "\033[0;34m[--xsecvar <cross-section var>]\033[0m                        \033[0;32mTThe variable to bin the cross section as a function of. Default is elec_E. Options are elec_E or elec_ang. \033[0m\n\n"
    "\033[0;34m[--xsecbins <cross-section bin mode> ]\033[0m                 \033[0;32mChoose whether to use a fine binned smearing/response matrix in truth. Default is standard. Options are standard or fine. \033[0m\n\n"
    "\033[0;34m[--xsec_smear <cross-section smear mode> ]\033[0m             \033[0;32mChoose how to define the cross section. Default is mcc8. Options are mcc8 or er. \033[0m\n\n"
    "-------------------------------------------------------\n\n";

    std::string usage3 = "\n\nTo run the detector systematics code, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --sys <systematics mode> [options (see below)]\033[0m \n\n"
    "\033[0;34m[--sys <systematics mode>]\033[0m                             \033[0;32mThe input mode of systematics to run. Options are default/ext/reweight \033[0m\n\n"
    "\033[0;34m[--xsec_smear <cross-section smear mode> ]\033[0m             \033[0;32mChoose how to define the cross section. Default is er. Options are mcc8 or er or wiener.\033[0m\n\n"
    "\033[0;34m[--binscaling <bin scaling option> ]\033[0m                   \033[0;32mChoose whether we want to apply a bin width scaling to the histogram. Default is width. Options are standard or width. \033[0m\n\n"
    "-------------------------------------------------------\n\n"
    "\n\nTo run the utility plotter code, run: \n\n"
    "\033[0;31m./nuexsec --run <run period num> --uplot <utility plotter mode>\033[0m \n\n"
    "\033[0;34m[--uplot <utility plotter mode>]\033[0m                        \033[0;32mThe input mode of the utility plotter to run. Options are default/bins/true/models. See UtilityPlotter Initalise function to see what functions are called. \033[0m\n\n"
    "This will run the utility plotting code\n\n";


    // -------------------------------------------------------------------------
    
    // Configure the utility class, this stores all the input configurations
    // that all classes will have access to
    _utility.Initalise(argc, argv, usage, usage2, usage3);

    // -------------------------------------------------------------------------

    // Initialise the selction script
    if (_utility.run_selection)  _selection_instance.Initialise(_utility);

    // Print the selection results
    if (_utility.print)          _phelper.Initialise(_utility );

    // Run the make histogram function
    if (_utility.make_histos)    _hplot.MakeHistograms(_utility);

    // Run the calculate cross section function
    if (_utility.calc_cross_sec) _xsec.Initialise(_utility);

    // Run the systematics helper code
    if (_utility.run_sys)        _syshelper.Initialise(_utility);

    // Run the utility plotting code
    if (_utility.run_uplot)      _uplot.Initialise(_utility);

    // -------------------------------------------------------------------------
    // Finished!
    std::cout << "\033[0;32m*** \t Exiting C++ Code... \t *** \033[0m" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();  // end time of script
    auto duration_sec = std::chrono::duration_cast<std::chrono::seconds>(stop - start); // time taken to run script
    auto duration_min = std::chrono::duration_cast<std::chrono::minutes>(stop - start); // time taken to run script
    std::cout << "Time taken by function: " << duration_sec.count() << " seconds" << std::endl; 
    std::cout << "Time taken by function: " << duration_min.count() << " minutes" << std::endl; 
    
    if (_utility.use_gpvm)
        std::cout << _utility.red << "If your terminal is not returned to you after 20 seconds of seeing this message then hit crtl c" << _utility.reset << std::endl; 
    
    // exit(0);
    return 0;
}
