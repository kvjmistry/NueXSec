#include "../include/PrintHelper.h"

// -----------------------------------------------------------------------------
void PrintHelper::Initialise(const char* run_period, const char * mc_file_in, bool _print_mc, bool _print_data, bool _print_ext, bool _print_dirt, Utility _utility ){

    std::cout << "Initalising Print Helper..." << std::endl;

    _util = _utility;

    std::string file_out_str = mc_file_in;

    std::string file_name;

    print_mc    = _print_mc;
    print_data  = _print_data;
    print_ext   = _print_ext;
    print_dirt  = _print_dirt;

    // define this to scale the MC to a desired POT
    // double additional_scaling = 6.0e20/_util.config_v.at(_util.k_Run1_Data_POT); // Use this to scale the POT to set amount -- also chnage the run POT
    double additional_scaling = 1.0; // Use this to use default
    if ( additional_scaling != 1.0) std::cout << "\033[0;34mWarning using an additional POT scale factor to print the selection results\033[0m" << std::endl;

    // Set the scale factors
    if (strcmp(run_period, "1") == 0){
        mc_scale_factor     = additional_scaling * _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_MC_POT);
        dirt_scale_factor   = additional_scaling * _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_Dirt_POT);
        intime_scale_factor = additional_scaling * _util.config_v.at(_util.k_Run1_Data_trig) / _util.config_v.at(_util.k_Run1_EXT_trig);
    }
    else if (strcmp(run_period, "3") == 0){
        mc_scale_factor     = additional_scaling * _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_MC_POT);
        dirt_scale_factor   = additional_scaling * _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_Dirt_POT);
        intime_scale_factor = additional_scaling * _util.config_v.at(_util.k_Run3_Data_trig) / _util.config_v.at(_util.k_Run3_EXT_trig);
    }
    else {
        std::cout << "Error Krish... You havent defined the run X POT numbers yet you donut!" << std::endl;
        exit(1);
    }
    
    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << "Scale Factors:\n" <<
    "MC Scale factor:        " << mc_scale_factor          << "\n" <<
    "Dirt Scale factor:      " << dirt_scale_factor        << "\n" <<
    "EXT Scale factor:       " << intime_scale_factor
    << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    // Data -----------------------
    if (print_data){
        file_name = Form("files/trees/nuexsec_selected_tree_data_run%s.root", run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_data = new TFile(file_name.c_str(), "READ");
            _util.GetTree(f_data, data_counter_tree, "data_counter_tree");
            tree_total_entries =  data_counter_tree->GetEntries();
        }
    }

    // EXT -----------------------
    if (print_ext){
        file_name = Form("files/trees/nuexsec_selected_tree_ext_run%s.root", run_period);
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_ext = new TFile(file_name.c_str(), "READ");
            _util.GetTree(f_ext, ext_counter_tree, "ext_counter_tree");
            tree_total_entries =  ext_counter_tree->GetEntries();
        }
    }

    // dirt -----------------------
    if (print_dirt){
        file_name = Form("files/trees/nuexsec_selected_tree_dirt_run%s.root", run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_dirt = new TFile(file_name.c_str(), "READ");
            _util.GetTree(f_dirt, dirt_counter_tree, "dirt_counter_tree");
            tree_total_entries =  dirt_counter_tree->GetEntries();
        }
    }

    // MC -----------------------
    if (print_mc){
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_mc_run%s.root", run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
            f_mc = new TFile( file_name.c_str(), "UPDATE");
            
            _util.GetTree(f_mc, mc_counter_tree, "mc_counter_tree");
            tree_total_entries =  mc_counter_tree->GetEntries();

            eff_tree  = new TTree("mc_eff_tree",     "mc_eff_tree");
            eff_tree->Branch("efficiency", &efficiency, "efficiency/D");
            eff_tree->Branch("purity", &purity, "purity/D");
        }
    }
    
    // Counter Tree
    if (print_mc){
        mc_counter_tree->SetBranchAddress("count_nue_cc_qe",  &count_nue_cc_qe);
        mc_counter_tree->SetBranchAddress("count_nue_cc_res", &count_nue_cc_res);
        mc_counter_tree->SetBranchAddress("count_nue_cc_dis", &count_nue_cc_dis);
        mc_counter_tree->SetBranchAddress("count_nue_cc_coh", &count_nue_cc_coh);
        mc_counter_tree->SetBranchAddress("count_nue_cc_mec", &count_nue_cc_mec);
        
        mc_counter_tree->SetBranchAddress("count_nuebar_cc_qe",  &count_nuebar_cc_qe);
        mc_counter_tree->SetBranchAddress("count_nuebar_cc_res", &count_nuebar_cc_res);
        mc_counter_tree->SetBranchAddress("count_nuebar_cc_dis", &count_nuebar_cc_dis);
        mc_counter_tree->SetBranchAddress("count_nuebar_cc_coh", &count_nuebar_cc_coh);
        mc_counter_tree->SetBranchAddress("count_nuebar_cc_mec", &count_nuebar_cc_mec);
        
        mc_counter_tree->SetBranchAddress("count_nue_cc_infv",      &count_nue_cc_infv);
        mc_counter_tree->SetBranchAddress("count_nuebar_cc_infv",   &count_nuebar_cc_infv);
        mc_counter_tree->SetBranchAddress("count_nue_cc_incryo",    &count_nue_cc_incryo);
        mc_counter_tree->SetBranchAddress("count_nuebar_cc_incryo", &count_nuebar_cc_incryo);
        
        mc_counter_tree->SetBranchAddress("count_numu_cc_qe",  &count_numu_cc_qe);
        mc_counter_tree->SetBranchAddress("count_numu_cc_res", &count_numu_cc_res);
        mc_counter_tree->SetBranchAddress("count_numu_cc_dis", &count_numu_cc_dis);
        mc_counter_tree->SetBranchAddress("count_numu_cc_coh", &count_numu_cc_coh);
        mc_counter_tree->SetBranchAddress("count_numu_cc_mec", &count_numu_cc_mec);
        
        mc_counter_tree->SetBranchAddress("count_numubar_cc_qe",  &count_numubar_cc_qe);
        mc_counter_tree->SetBranchAddress("count_numubar_cc_res", &count_numubar_cc_res);
        mc_counter_tree->SetBranchAddress("count_numubar_cc_dis", &count_numubar_cc_dis);
        mc_counter_tree->SetBranchAddress("count_numubar_cc_coh", &count_numubar_cc_coh);
        mc_counter_tree->SetBranchAddress("count_numubar_cc_mec", &count_numubar_cc_mec);
        
        mc_counter_tree->SetBranchAddress("count_numu_cc_infv",      &count_numu_cc_infv);
        mc_counter_tree->SetBranchAddress("count_numubar_cc_infv",   &count_numubar_cc_infv);
        mc_counter_tree->SetBranchAddress("count_numu_cc_incryo",    &count_numu_cc_incryo);
        mc_counter_tree->SetBranchAddress("count_numubar_cc_incryo", &count_numubar_cc_incryo);

        mc_counter_tree->SetBranchAddress("count_pi0_nue_cc_nopi0",    &count_pi0_nue_cc_nopi0);
        mc_counter_tree->SetBranchAddress("count_pi0_nue_cc_pi0",      &count_pi0_nue_cc_pi0);
        mc_counter_tree->SetBranchAddress("count_pi0_nuebar_cc_nopi0", &count_pi0_nuebar_cc_nopi0);
        mc_counter_tree->SetBranchAddress("count_pi0_nuebar_cc_pi0",   &count_pi0_nuebar_cc_pi0);
        mc_counter_tree->SetBranchAddress("count_pi0_numu_cc_nopi0",   &count_pi0_numu_cc_nopi0);
        mc_counter_tree->SetBranchAddress("count_pi0_numu_cc_pi0",     &count_pi0_numu_cc_pi0);
        mc_counter_tree->SetBranchAddress("count_pi0_nc_nopi0",        &count_pi0_nc_nopi0);
        mc_counter_tree->SetBranchAddress("count_pi0_nc_pi0",          &count_pi0_nc_pi0);
        
        mc_counter_tree->SetBranchAddress("count_nue_cc",          &count_nue_cc);
        mc_counter_tree->SetBranchAddress("count_nuebar_cc",       &count_nuebar_cc);
        mc_counter_tree->SetBranchAddress("count_nu_out_fv",       &count_nu_out_fv);
        mc_counter_tree->SetBranchAddress("count_cosmic",          &count_cosmic);
        mc_counter_tree->SetBranchAddress("count_numu_cc",         &count_numu_cc);
        mc_counter_tree->SetBranchAddress("count_numu_cc_pi0",     &count_numu_cc_pi0);
        mc_counter_tree->SetBranchAddress("count_nc",              &count_nc);
        mc_counter_tree->SetBranchAddress("count_nc_pi0",          &count_nc_pi0);
        mc_counter_tree->SetBranchAddress("count_unmatched",       &count_unmatched);
        mc_counter_tree->SetBranchAddress("count_unmatched_nue",   &count_unmatched_nue);
        mc_counter_tree->SetBranchAddress("count_cosmic_nue",      &count_cosmic_nue);
        mc_counter_tree->SetBranchAddress("count_unmatched_nuebar",&count_unmatched_nuebar);
        mc_counter_tree->SetBranchAddress("count_cosmic_nuebar",   &count_cosmic_nuebar);
        mc_counter_tree->SetBranchAddress("count_total_mc",        &count_total_mc);
    }

    if (print_data){
        data_counter_tree->SetBranchAddress("count_data",         &count_data);
    }
    
    if (print_ext){
        ext_counter_tree->SetBranchAddress("count_ext",          &count_ext);
    }
    
    if (print_dirt){
        dirt_counter_tree->SetBranchAddress("count_dirt",         &count_dirt);
    }

    PrintResults();

}
// -----------------------------------------------------------------------------
void PrintHelper::PrintResults(){

    // Loop over the cuts
    for (int p =0; p < tree_total_entries; p++){

        if (print_mc)   mc_counter_tree->GetEntry(p);
        if (print_data) data_counter_tree->GetEntry(p);
        if (print_ext)  ext_counter_tree->GetEntry(p);
        if (print_dirt) dirt_counter_tree->GetEntry(p);

        // Set counters at the start
        if (p == 0) {
            tot_true_infv_nues = count_nue_cc + count_nuebar_cc + count_unmatched_nue + count_cosmic_nue + count_unmatched_nuebar + count_cosmic_nuebar;

            init_count_nue_cc       = count_nue_cc;
            init_count_nuebar_cc    = count_nuebar_cc;
            init_count_nu_out_fv    = count_nu_out_fv;
            init_count_cosmic       = count_cosmic;
            init_count_numu_cc      = count_numu_cc;
            init_count_numu_cc_pi0  = count_numu_cc_pi0;
            init_count_nc           = count_nc;
            init_count_nc_pi0       = count_nc_pi0;
            init_count_unmatched    = count_unmatched;
            init_count_unmatched_nue    = count_unmatched_nue;
            init_count_cosmic_nue       = count_cosmic_nue;
            init_count_unmatched_nuebar = count_unmatched_nuebar;
            init_count_cosmic_nuebar    = count_cosmic_nuebar;
            if (print_ext) init_count_ext = count_ext;
            if (print_dirt) init_count_dirt = count_dirt;

        }

        
        // Sum of selected mc, dirt and ext. The dirt and ext are scaled to the MC POT
        double sum_mc_dirt_ext = count_total_mc+ (count_ext * (intime_scale_factor / mc_scale_factor)) + (count_dirt * (dirt_scale_factor / mc_scale_factor));

        // Set the efficiency and purity for case zero
        if (p == 0){
            efficiency_last = double(count_nue_cc + count_nuebar_cc) / double(tot_true_infv_nues);
            purity_last     = double(count_nue_cc + count_nuebar_cc) / double(sum_mc_dirt_ext);
        }


        if (print_mc && print_data && print_ext && print_dirt){
            std::cout << "\n------------------------------------------------" << std::endl;
        }
        std::cout << "------------------------------------------------" << std::endl;
        
        std::cout.precision(6);
        std::cout << "\n\033[0;33m <" << _util.cut_dirs.at(p) << "> \033[0m" << std::endl;
        
        if (!print_mc) {
            if (p == 0) printf (" %-21s %-10s %-10s %-10s %-10s %-12s\n", " " , "MC POT", "Data POT", "Unscaled", "\% Remaining", " \% Change from Prev Cut");
        }
        else printf (" %-21s %-10s %-10s %-10s %-10s %-12s\n", " " , "MC POT", "Data POT", "Unscaled", "\% Remaining", " \% Change from Prev Cut");

        if (print_mc){
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Nue CC",               count_nue_cc,           double(count_nue_cc           * mc_scale_factor  ), " ", double( 100 * count_nue_cc / init_count_nue_cc),                     double(-100 * (count_nue_cc           - prev_count_nue_cc)           / prev_count_nue_cc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "NueBar CC",            count_nuebar_cc,        double(count_nuebar_cc        * mc_scale_factor  ), " ", double( 100 * count_nuebar_cc / init_count_nuebar_cc),               double(-100 * (count_nuebar_cc        - prev_count_nuebar_cc)        / prev_count_nuebar_cc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Nu out FV",            count_nu_out_fv,        double(count_nu_out_fv        * mc_scale_factor  ), " ", double( 100 * count_nu_out_fv / init_count_nu_out_fv),               double(-100 * (count_nu_out_fv        - prev_count_nu_out_fv)        / prev_count_nu_out_fv));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Cosmic",               count_cosmic,           double(count_cosmic           * mc_scale_factor  ), " ", double( 100 * count_cosmic / init_count_cosmic),                     double(-100 * (count_cosmic           - prev_count_cosmic)           / prev_count_cosmic));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Numu CC",              count_numu_cc,          double(count_numu_cc          * mc_scale_factor  ), " ", double( 100 * count_numu_cc / init_count_numu_cc),                   double(-100 * (count_numu_cc          - prev_count_numu_cc)          / prev_count_numu_cc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Numu CC Pi0",          count_numu_cc_pi0,      double(count_numu_cc_pi0      * mc_scale_factor  ), " ", double( 100 * count_numu_cc_pi0 / init_count_numu_cc_pi0),           double(-100 * (count_numu_cc_pi0      - prev_count_numu_cc_pi0)      / prev_count_numu_cc_pi0));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "NC",                   count_nc,               double(count_nc               * mc_scale_factor  ), " ", double( 100 * count_nc / init_count_nc),                             double(-100 * (count_nc               - prev_count_nc)               / prev_count_nc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "NC Pi0",               count_nc_pi0,           double(count_nc_pi0           * mc_scale_factor  ), " ", double( 100 * count_nc_pi0 / init_count_nc_pi0),                     double(-100 * (count_nc_pi0           - prev_count_nc_pi0)           / prev_count_nc_pi0));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Cosmic Nue CC",        count_cosmic_nue,       double(count_cosmic_nue       * mc_scale_factor  ), " ", double( 100 * count_cosmic_nue / init_count_cosmic_nue),             double(-100 * (count_cosmic_nue       - prev_count_cosmic_nue)       / prev_count_cosmic_nue));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Cosmic Nuebar CC",     count_cosmic_nuebar,    double(count_cosmic_nuebar    * mc_scale_factor  ), " ", double( 100 * count_cosmic_nuebar),                                  double(-100 * (count_cosmic_nuebar    - prev_count_cosmic_nuebar)    / prev_count_cosmic_nuebar));
            if (p == 0) printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Non-Reco'd Nue CC",    count_unmatched_nue,    double(count_unmatched_nue    * mc_scale_factor  ), " ", double( 100 * count_unmatched_nue / init_count_unmatched_nue),       double(-100 * (count_unmatched_nue    - prev_count_unmatched_nue)    / prev_count_unmatched_nue));
            if (p == 0) printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Non-Reco'd Nuebar CC", count_unmatched_nuebar, double(count_unmatched_nuebar * mc_scale_factor  ), " ", double( 100 * count_unmatched_nuebar / init_count_unmatched_nuebar), double(-100 * (count_unmatched_nuebar - prev_count_unmatched_nuebar) / prev_count_unmatched_nuebar));
       }

        if (print_dirt){
            printf (" %-20s: %-10.2f %-10.2f %-10.2f %f %9.2f\n", "Dirt", double(count_dirt * (dirt_scale_factor / mc_scale_factor)), double(count_dirt * dirt_scale_factor), count_dirt, double( 100 * count_dirt / init_count_dirt), double(100 * (prev_count_dirt - count_dirt) / prev_count_dirt) );
            
        }

        if (print_ext) { 
            printf (" %-20s: %-10.2f %-10.2f %-10.2f %f %9.2f\n", "Off-Beam Data", double(count_ext * (intime_scale_factor / mc_scale_factor)), double(count_ext * intime_scale_factor), count_ext, double( 100 * count_ext / init_count_ext), double( 100 * (prev_count_ext - count_ext) / prev_count_ext) );
            
        }

        if (print_mc && print_dirt){
            double tot_mc_bkg = count_numu_cc + count_numu_cc_pi0 + count_nc + count_nc_pi0;
            printf ("\n %-20s: %-10.2f %-10.2f\n", "Total Beam Bkg", tot_mc_bkg + count_dirt * (dirt_scale_factor / mc_scale_factor), double(tot_mc_bkg * mc_scale_factor + count_dirt * dirt_scale_factor ) );
        }

        if (print_mc && print_ext){
            double tot_cosmic_bkg = count_cosmic + count_cosmic_nue + count_cosmic_nuebar;
            printf ("\n %-20s: %-10.2f %-10.2f\n", "Total Cosmic Bkg", tot_cosmic_bkg + count_ext * (intime_scale_factor / mc_scale_factor), double(tot_cosmic_bkg * mc_scale_factor + count_ext * intime_scale_factor ));
        }
        
        if (print_mc){
            printf ("\n %-20s: %-10.2f %-10.2f\n", "Total Candidate Nue", sum_mc_dirt_ext, double(sum_mc_dirt_ext         * mc_scale_factor  ));
        }
  
        if (print_mc){
            std::cout << "\n----------- Neutrinos in FV Truth -------------" << std::endl;
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Nue CC QE",    count_nue_cc_qe*mc_scale_factor,  "Nuebar CC QE",    count_nuebar_cc_qe*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Nue CC Res",   count_nue_cc_res*mc_scale_factor, "Nuebar CC Res",   count_nuebar_cc_res*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Nue CC DIS",   count_nue_cc_dis*mc_scale_factor, "Nuebar CC DIS",   count_nuebar_cc_dis*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Nue CC COH",   count_nue_cc_coh*mc_scale_factor, "Nuebar CC COH",   count_nuebar_cc_coh*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Nue CC MEC",   count_nue_cc_mec*mc_scale_factor, "Nuebar CC MEC",   count_nuebar_cc_mec);
            std::cout << std::endl;
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Numu CC QE",    count_numu_cc_qe*mc_scale_factor,  "Numubar CC QE",    count_numubar_cc_qe*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Numu CC Res",   count_numu_cc_res*mc_scale_factor, "Numubar CC Res",   count_numubar_cc_res*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Numu CC DIS",   count_numu_cc_dis*mc_scale_factor, "Numubar CC DIS",   count_numubar_cc_dis*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Numu CC COH",   count_numu_cc_coh*mc_scale_factor, "Numubar CC COH",   count_numubar_cc_coh*mc_scale_factor);
            printf (" %-12s: %-10.2f %-12s: %-10.2f\n", "Numu CC MEC",   count_numu_cc_mec*mc_scale_factor, "Numubar CC MEC",   count_numubar_cc_mec*mc_scale_factor);
            std::cout << "------------------------------------------------" << std::endl;
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot Nue in FV",     count_nue_cc_infv*mc_scale_factor,     "Tot Nue in Cryo",       count_nue_cc_incryo*mc_scale_factor);
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot Nuebar in FV",  count_nuebar_cc_infv*mc_scale_factor,  "Tot Nuebar in Cryo",    count_nuebar_cc_incryo*mc_scale_factor);
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot NuMu in FV",    count_numu_cc_infv*mc_scale_factor,    "Tot NuMu in Cryo",      count_numu_cc_incryo*mc_scale_factor);
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot NuMubar in FV", count_numubar_cc_infv*mc_scale_factor, "Tot NuMubar in Cryo",   count_numubar_cc_incryo*mc_scale_factor);
            std::cout << "------------------------------------------------" << std::endl;
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot Nue 0 Pi0",     count_pi0_nue_cc_nopi0*mc_scale_factor,     "Tot Nuebar 0 Pi0",  count_pi0_nuebar_cc_nopi0*mc_scale_factor);
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot Nue Pi0",       count_pi0_nue_cc_pi0*mc_scale_factor,       "Tot Nuebar Pi0",    count_pi0_nuebar_cc_pi0*mc_scale_factor);
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot Numu 0 Pi0",    count_pi0_numu_cc_nopi0*mc_scale_factor,    "Tot NC 0 Pi0",      count_pi0_nc_nopi0*mc_scale_factor);
            printf (" %-17s: %-10.2f %-20s: %-10.2f\n", "Tot Numu Pi0",      count_pi0_numu_cc_pi0*mc_scale_factor,      "Tot NC Pi0",        count_pi0_nc_pi0*mc_scale_factor);
        }
        
        if (print_mc){
            std::cout << "------------------------------------------------" << std::endl;
            efficiency = double(count_nue_cc + count_nuebar_cc) / double(tot_true_infv_nues);
            purity     = double(count_nue_cc + count_nuebar_cc) / double(sum_mc_dirt_ext);
            printf (" %-15s: ( %-6.1f / %-7.1f ) = %-3.2f %% \n", "Purity", count_nue_cc + count_nuebar_cc, sum_mc_dirt_ext, 100 * purity);
            printf (" %-15s: ( %-6.1f / %-7.1f ) = %-3.2f %% \n", "Efficiency", count_nue_cc + count_nuebar_cc, tot_true_infv_nues, 100 * efficiency);
            
            std::cout << std::endl;
            std::cout << " Purity Change     : " <<  100 * (purity - purity_last) << " \%" << std::endl;
            std::cout << " Efficiency Change : "  << 100 * (efficiency - efficiency_last) << " \%" << std::endl;

            efficiency_last = efficiency;
            purity_last     = purity;
        }

        if (print_data){
            std::cout << "------------------------------------------------" << std::endl;
            double tot_bkg = 0;
            if (print_mc && print_dirt && print_ext) tot_bkg = (count_nu_out_fv + count_cosmic + count_numu_cc + count_numu_cc_pi0 + count_nc + count_nc_pi0 + count_unmatched) * mc_scale_factor +
                              count_dirt * dirt_scale_factor + count_ext * intime_scale_factor;
            
            std::cout << " Total Selected Data : " << count_data << std::endl;
            
            if (print_mc && print_dirt && print_ext) std::cout << " Total Nue Candidates in data : " << count_data - tot_bkg<< std::endl;
            
            std::cout << "------------------------------------------------" << std::endl;
        }

        // Fill the efficiency tree
        if (print_mc && print_data && print_ext && print_dirt) eff_tree->Fill();

        // Set counters for the previous cut
        prev_count_nue_cc       = count_nue_cc;
        prev_count_nuebar_cc    = count_nuebar_cc;
        prev_count_nu_out_fv    = count_nu_out_fv;
        prev_count_cosmic       = count_cosmic;
        prev_count_numu_cc      = count_numu_cc;
        prev_count_numu_cc_pi0  = count_numu_cc_pi0;
        prev_count_nc           = count_nc;
        prev_count_nc_pi0       = count_nc_pi0;
        prev_count_unmatched    = count_unmatched;
        prev_count_unmatched_nue    = count_unmatched_nue;
        prev_count_cosmic_nue       = count_cosmic_nue;
        prev_count_unmatched_nuebar = count_unmatched_nuebar;
        prev_count_cosmic_nuebar    = count_cosmic_nuebar;
        if (print_ext) prev_count_ext = count_ext;
        if (print_dirt) prev_count_dirt = count_dirt;
        
    }

    // Save the efficiency tree
    if (print_mc && print_data && print_ext && print_dirt) {
        
        if (print_mc){
            f_mc->cd();
            eff_tree->Write("",TObject::kOverwrite);
        }
        
    }

}
