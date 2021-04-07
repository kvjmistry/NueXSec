#include "../include/PrintHelper.h"

// -----------------------------------------------------------------------------
void PrintHelper::Initialise(Utility _utility ){

    std::cout << "Initalising Print Helper..." << std::endl;

    _util = _utility;

    std::string file_out_str = std::string(_util.mc_file_name);

    std::string file_name;

    // define this to scale the MC to a desired POT
    // double additional_scaling = 2.0e20/_util.config_v.at(_util.k_Run1_Data_POT); // Use this to scale the POT to set amount -- also chnage the run POT
    double additional_scaling = 1.0; // Use this to use default
    if ( additional_scaling != 1.0) std::cout << "\033[0;34mWarning using an additional POT scale factor to print the selection results\033[0m" << std::endl;

    // Set the scale factors
    mc_scale_factor     = additional_scaling * _util.mc_scale_factor;
    dirt_scale_factor   = additional_scaling * _util.dirt_scale_factor;
    ext_scale_factor    = additional_scaling * _util.ext_scale_factor;
    
    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << "Scale Factors:\n" <<
    "MC Scale factor:        " << mc_scale_factor          << "\n" <<
    "Dirt Scale factor:      " << dirt_scale_factor        << "\n" <<
    "EXT Scale factor:       " << ext_scale_factor
    << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    // Data -----------------------
    if (_util.print_data){
        file_name = Form("files/trees/nuexsec_selected_tree_data_run%s.root", _util.run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_data = TFile::Open(file_name.c_str(), "READ");
            _util.GetTree(f_data, data_counter_tree, "data_counter_tree");
            tree_total_entries =  data_counter_tree->GetEntries();
        }

        // Get the histogram files too
        file_name = Form("files/nuexsec_data_run%s.root", _util.run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_data_hist = TFile::Open(file_name.c_str(), "READ");
        }

    }

    // EXT -----------------------
    if (_util.print_ext){
        file_name = Form("files/trees/nuexsec_selected_tree_ext_run%s.root", _util.run_period);
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_ext = TFile::Open(file_name.c_str(), "READ");
            _util.GetTree(f_ext, ext_counter_tree, "ext_counter_tree");
            tree_total_entries =  ext_counter_tree->GetEntries();
        }

        // Get the histogram files too
        file_name = Form("files/nuexsec_ext_run%s.root", _util.run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_ext_hist = TFile::Open(file_name.c_str(), "READ");
        }
    }

    // dirt -----------------------
    if (_util.print_dirt){
        file_name = Form("files/trees/nuexsec_selected_tree_dirt_run%s.root", _util.run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_dirt = TFile::Open(file_name.c_str(), "READ");
            _util.GetTree(f_dirt, dirt_counter_tree, "dirt_counter_tree");
            tree_total_entries =  dirt_counter_tree->GetEntries();
        }

        // Get the histogram files too
        file_name = Form("files/nuexsec_dirt_run%s.root", _util.run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_dirt_hist = TFile::Open(file_name.c_str(), "READ");
        }
    }

    // MC -----------------------
    if (_util.print_mc){
        if (file_out_str == "empty") file_name = Form("files/trees/nuexsec_selected_tree_mc_run%s.root", _util.run_period);
        else file_name = "files/trees/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
            f_mc = TFile::Open( file_name.c_str(), "UPDATE");
            
            _util.GetTree(f_mc, mc_counter_tree, "mc_counter_tree");
            _util.GetTree(f_mc, mc_nue_counter_tree, "mc_nue_counter_tree");
            tree_total_entries =  mc_counter_tree->GetEntries();

            eff_tree  = new TTree("mc_eff_tree",     "mc_eff_tree");
            eff_tree->SetDirectory(f_mc);
            eff_tree->Branch("efficiency", &efficiency, "efficiency/D");
            eff_tree->Branch("eff_err", &eff_err, "eff_err/D");
            eff_tree->Branch("purity", &purity, "purity/D");
        }

        // Get the histogram files too
        file_name = Form("files/nuexsec_mc_run%s.root", _util.run_period);

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_mc_hist = TFile::Open(file_name.c_str(), "READ");
        }
    }
    
    // Counter Tree
    if (_util.print_mc){
        mc_nue_counter_tree->SetBranchAddress("count_nue_cc_qe",  &count_nue_cc_qe);
        mc_nue_counter_tree->SetBranchAddress("count_nue_cc_res", &count_nue_cc_res);
        mc_nue_counter_tree->SetBranchAddress("count_nue_cc_dis", &count_nue_cc_dis);
        mc_nue_counter_tree->SetBranchAddress("count_nue_cc_coh", &count_nue_cc_coh);
        mc_nue_counter_tree->SetBranchAddress("count_nue_cc_mec", &count_nue_cc_mec);
        
        mc_nue_counter_tree->SetBranchAddress("count_nuebar_cc_qe",  &count_nuebar_cc_qe);
        mc_nue_counter_tree->SetBranchAddress("count_nuebar_cc_res", &count_nuebar_cc_res);
        mc_nue_counter_tree->SetBranchAddress("count_nuebar_cc_dis", &count_nuebar_cc_dis);
        mc_nue_counter_tree->SetBranchAddress("count_nuebar_cc_coh", &count_nuebar_cc_coh);
        mc_nue_counter_tree->SetBranchAddress("count_nuebar_cc_mec", &count_nuebar_cc_mec);
        
        mc_nue_counter_tree->SetBranchAddress("count_nue_cc_infv",      &count_nue_cc_infv);
        mc_nue_counter_tree->SetBranchAddress("count_nuebar_cc_infv",   &count_nuebar_cc_infv);
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

        mc_nue_counter_tree->SetBranchAddress("count_pi0_nue_cc_nopi0",    &count_pi0_nue_cc_nopi0);
        mc_nue_counter_tree->SetBranchAddress("count_pi0_nue_cc_pi0",      &count_pi0_nue_cc_pi0);
        mc_nue_counter_tree->SetBranchAddress("count_pi0_nuebar_cc_nopi0", &count_pi0_nuebar_cc_nopi0);
        mc_nue_counter_tree->SetBranchAddress("count_pi0_nuebar_cc_pi0",   &count_pi0_nuebar_cc_pi0);
        mc_counter_tree    ->SetBranchAddress("count_pi0_numu_cc_nopi0",   &count_pi0_numu_cc_nopi0);
        mc_counter_tree    ->SetBranchAddress("count_pi0_numu_cc_pi0",     &count_pi0_numu_cc_pi0);
        mc_counter_tree    ->SetBranchAddress("count_pi0_nc_nopi0",        &count_pi0_nc_nopi0);
        mc_counter_tree    ->SetBranchAddress("count_pi0_nc_pi0",          &count_pi0_nc_pi0);
        
        mc_nue_counter_tree->SetBranchAddress("count_nue_cc",          &count_nue_cc);
        mc_nue_counter_tree->SetBranchAddress("count_nuebar_cc",       &count_nuebar_cc);
        mc_counter_tree    ->SetBranchAddress("count_nu_out_fv",       &count_nu_out_fv);
        mc_counter_tree    ->SetBranchAddress("count_cosmic",          &count_cosmic);
        mc_counter_tree    ->SetBranchAddress("count_numu_cc",         &count_numu_cc);
        mc_counter_tree    ->SetBranchAddress("count_numu_cc_pi0",     &count_numu_cc_pi0);
        mc_counter_tree    ->SetBranchAddress("count_nc",              &count_nc);
        mc_counter_tree    ->SetBranchAddress("count_nc_pi0",          &count_nc_pi0);
        mc_counter_tree    ->SetBranchAddress("count_unmatched",       &count_unmatched);
        mc_nue_counter_tree->SetBranchAddress("count_unmatched_nue",   &count_unmatched_nue);
        mc_nue_counter_tree->SetBranchAddress("count_cosmic_nue",      &count_cosmic_nue);
        mc_nue_counter_tree->SetBranchAddress("count_unmatched_nuebar",&count_unmatched_nuebar);
        mc_nue_counter_tree->SetBranchAddress("count_cosmic_nuebar",   &count_cosmic_nuebar);
        mc_nue_counter_tree->SetBranchAddress("count_thr_nue",         &count_thr_nue);
        mc_nue_counter_tree->SetBranchAddress("count_thr_nuebar",      &count_thr_nuebar);
        mc_counter_tree    ->SetBranchAddress("count_total_mc",        &count_total_mc);
    }

    if (_util.print_data){
        data_counter_tree->SetBranchAddress("count_data",         &count_data);
    }
    
    if (_util.print_ext){
        ext_counter_tree->SetBranchAddress("count_ext",          &count_ext);
    }
    
    if (_util.print_dirt){
        dirt_counter_tree->SetBranchAddress("count_dirt",         &count_dirt);
    }

    // Lets get the histograms here
    if (_util.print_mc) GetHists();


    PrintResults();

}
// -----------------------------------------------------------------------------
void PrintHelper::PrintResults(){

    // Loop over the cuts
    for (int p =0; p < tree_total_entries; p++){

        if (_util.print_mc)   mc_counter_tree->GetEntry(p);
        if (_util.print_mc)   mc_nue_counter_tree->GetEntry(p);
        if (_util.print_data) data_counter_tree->GetEntry(p);
        if (_util.print_ext)  ext_counter_tree->GetEntry(p);
        if (_util.print_dirt) dirt_counter_tree->GetEntry(p);

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
            init_count_thr_nue          = count_thr_nue;
            init_count_thr_nuebar       = count_thr_nuebar;
            if (_util.print_ext) init_count_ext = count_ext;
            if (_util.print_dirt) init_count_dirt = count_dirt;

            init_count_nue_cc_qe  = count_nue_cc_qe;
            init_count_nue_cc_res = count_nue_cc_res;
            init_count_nue_cc_dis = count_nue_cc_dis;
            init_count_nue_cc_coh = count_nue_cc_coh;
            init_count_nue_cc_mec = count_nue_cc_mec;

            init_count_nuebar_cc_qe  = count_nuebar_cc_qe;
            init_count_nuebar_cc_res = count_nuebar_cc_res;
            init_count_nuebar_cc_dis = count_nuebar_cc_dis;
            init_count_nuebar_cc_coh = count_nuebar_cc_coh;
            init_count_nuebar_cc_mec = count_nuebar_cc_mec;

        }

        
        // Sum of selected mc, dirt and ext. The dirt and ext are scaled to the MC POT
        double sum_mc_dirt_ext = count_total_mc+ (count_ext * (ext_scale_factor / mc_scale_factor)) + (count_dirt * (dirt_scale_factor / mc_scale_factor));

        // Set the efficiency and purity for case zero
        if (p == 0){
            efficiency_last = double(count_nue_cc + count_nuebar_cc) / double(tot_true_infv_nues);
            purity_last     = double(count_nue_cc + count_nuebar_cc) / double(sum_mc_dirt_ext);
        }


        if (_util.print_mc && _util.print_data && _util.print_ext && _util.print_dirt){
            std::cout << "\n------------------------------------------------" << std::endl;
        }
        std::cout << "------------------------------------------------" << std::endl;
        
        std::cout.precision(6);
        std::cout << "\n\033[0;33m <" << _util.cut_dirs.at(p) << "> \033[0m" << std::endl;
        
        if (!_util.print_mc) {
            if (p == 0) printf (" %-21s %-10s %-10s %-10s %-10s %-12s\n", " " , "MC POT", "Data POT", "Unscaled", "\% Remaining", " \% Change from Prev Cut");
        }
        else printf (" %-21s %-10s %-10s %-10s %-10s %-12s\n", " " , "MC POT", "Data POT", "Unscaled", "\% Remaining", " \% Change from Prev Cut");

        if (_util.print_mc){
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Nue CC",               count_nue_cc,           double(count_nue_cc           * mc_scale_factor  ), " ", double( 100 * count_nue_cc / init_count_nue_cc),                     double(-100 * (count_nue_cc           - prev_count_nue_cc)           / prev_count_nue_cc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "NueBar CC",            count_nuebar_cc,        double(count_nuebar_cc        * mc_scale_factor  ), " ", double( 100 * count_nuebar_cc / init_count_nuebar_cc),               double(-100 * (count_nuebar_cc        - prev_count_nuebar_cc)        / prev_count_nuebar_cc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Nu out FV",            count_nu_out_fv,        double(count_nu_out_fv        * mc_scale_factor  ), " ", double( 100 * count_nu_out_fv / init_count_nu_out_fv),               double(-100 * (count_nu_out_fv        - prev_count_nu_out_fv)        / prev_count_nu_out_fv));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Numu CC",              count_numu_cc,          double(count_numu_cc          * mc_scale_factor  ), " ", double( 100 * count_numu_cc / init_count_numu_cc),                   double(-100 * (count_numu_cc          - prev_count_numu_cc)          / prev_count_numu_cc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Numu CC Pi0",          count_numu_cc_pi0,      double(count_numu_cc_pi0      * mc_scale_factor  ), " ", double( 100 * count_numu_cc_pi0 / init_count_numu_cc_pi0),           double(-100 * (count_numu_cc_pi0      - prev_count_numu_cc_pi0)      / prev_count_numu_cc_pi0));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "NC",                   count_nc,               double(count_nc               * mc_scale_factor  ), " ", double( 100 * count_nc / init_count_nc),                             double(-100 * (count_nc               - prev_count_nc)               / prev_count_nc));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "NC Pi0",               count_nc_pi0,           double(count_nc_pi0           * mc_scale_factor  ), " ", double( 100 * count_nc_pi0 / init_count_nc_pi0),                     double(-100 * (count_nc_pi0           - prev_count_nc_pi0)           / prev_count_nc_pi0));
            
            if (_util.print_dirt){
                printf (" %-20s: %-10.2f %-10.2f %-10.2f %f %9.2f\n", "Dirt", double(count_dirt * (dirt_scale_factor / mc_scale_factor)), double(count_dirt * dirt_scale_factor), count_dirt, double( 100 * count_dirt / init_count_dirt), double(100 * (prev_count_dirt - count_dirt) / prev_count_dirt) );
            }
            
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Cosmic",               count_cosmic,           double(count_cosmic           * mc_scale_factor  ), " ", double( 100 * count_cosmic / init_count_cosmic),                     double(-100 * (count_cosmic           - prev_count_cosmic)           / prev_count_cosmic));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Cosmic Nue CC",        count_cosmic_nue,       double(count_cosmic_nue       * mc_scale_factor  ), " ", double( 100 * count_cosmic_nue / init_count_cosmic_nue),             double(-100 * (count_cosmic_nue       - prev_count_cosmic_nue)       / prev_count_cosmic_nue));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Cosmic Nuebar CC",     count_cosmic_nuebar,    double(count_cosmic_nuebar    * mc_scale_factor  ), " ", double( 100 * count_cosmic_nuebar),                                  double(-100 * (count_cosmic_nuebar    - prev_count_cosmic_nuebar)    / prev_count_cosmic_nuebar));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Below Th Nue CC",      count_thr_nue,          double(count_thr_nue          * mc_scale_factor  ), " ", double( 100 * count_thr_nue),                                        double(-100 * (count_thr_nue       - prev_count_thr_nue)       / prev_count_thr_nue));
            printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Below Th Nuebar CC",   count_thr_nuebar,       double(count_thr_nuebar       * mc_scale_factor  ), " ", double( 100 * count_thr_nuebar),                                     double(-100 * (count_thr_nuebar    - prev_count_thr_nuebar)    / prev_count_thr_nuebar));
            if (p == 0) printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Non-Reco'd Nue CC",    count_unmatched_nue,    double(count_unmatched_nue    * mc_scale_factor  ), " ", double( 100 * count_unmatched_nue / init_count_unmatched_nue),       double(-100 * (count_unmatched_nue    - prev_count_unmatched_nue)    / prev_count_unmatched_nue));
            if (p == 0) printf (" %-20s: %-10.2f %-10.2f %-10s %-12.1f %-10.1f\n", "Non-Reco'd Nuebar CC", count_unmatched_nuebar, double(count_unmatched_nuebar * mc_scale_factor  ), " ", double( 100 * count_unmatched_nuebar / init_count_unmatched_nuebar), double(-100 * (count_unmatched_nuebar - prev_count_unmatched_nuebar) / prev_count_unmatched_nuebar));
       }

        if (_util.print_ext) { 
            printf (" %-20s: %-10.2f %-10.2f %-10.2f %f %9.2f\n", "Off-Beam Data", double(count_ext * (ext_scale_factor / mc_scale_factor)), double(count_ext * ext_scale_factor), count_ext, double( 100 * count_ext / init_count_ext), double( 100 * (prev_count_ext - count_ext) / prev_count_ext) );
            
        }

        if (_util.print_mc && _util.print_dirt){
            double tot_mc_bkg = count_nu_out_fv + count_numu_cc + count_numu_cc_pi0 + count_nc + count_nc_pi0 + count_thr_nue + count_thr_nuebar + count_cosmic + count_cosmic_nue + count_cosmic_nuebar;
            // printf ("\n %-20s: %-10.2f %-10.2f\n", "Total Beam Bkg", tot_mc_bkg + count_dirt * (dirt_scale_factor / mc_scale_factor), double(tot_mc_bkg * mc_scale_factor + count_dirt * dirt_scale_factor ) );
            printf ("\n %-20s: %-10.2f %-10.2f\n", "Total Overlay Bkg", tot_mc_bkg, double(tot_mc_bkg * mc_scale_factor ) );
        }

        if (_util.print_mc && _util.print_ext){
            double tot_cosmic_bkg = count_cosmic + count_cosmic_nue + count_cosmic_nuebar;
            printf ("\n %-20s: %-10.2f %-10.2f\n", "Total Cosmic Bkg", tot_cosmic_bkg + count_ext * (ext_scale_factor / mc_scale_factor), double(tot_cosmic_bkg * mc_scale_factor + count_ext * ext_scale_factor ));
        }

        if (_util.print_mc){
            double tot_MC = count_nue_cc + count_nuebar_cc + count_nu_out_fv + count_cosmic + count_numu_cc + count_numu_cc_pi0 + count_nc + count_nc_pi0 + count_cosmic_nue + count_cosmic_nuebar + count_thr_nue + count_thr_nuebar;
            printf ("\n %-20s: %-10.2f %-10.2f\n", "Total Overlay S+B", tot_MC , double(tot_MC         * mc_scale_factor  ));
        }
        
        if (_util.print_mc){
            printf ("\n %-20s: %-10.2f %-10.2f\n", "Total MC+EXT+Dirt", sum_mc_dirt_ext, double(sum_mc_dirt_ext         * mc_scale_factor  ));
        }
  
        if (_util.print_mc){
            std::cout << "\n----------- Neutrinos in FV Truth -------------" << std::endl;
            printf (" %-12s: %-7.2f(%2.1f%%) %-12s: %-8.2f(%2.1f%%)\n", "Nue CC QE ",   count_nue_cc_qe*mc_scale_factor , 100 * count_nue_cc_qe  / init_count_nue_cc_qe,  "Nuebar CC QE",    count_nuebar_cc_qe*mc_scale_factor,  100 * count_nuebar_cc_qe  / init_count_nuebar_cc_qe);
            printf (" %-12s: %-7.2f(%2.1f%%) %-12s: %-7.2f(%2.1f%%)\n", "Nue CC Res",   count_nue_cc_res*mc_scale_factor, 100 * count_nue_cc_res / init_count_nue_cc_res, "Nuebar CC Res",   count_nuebar_cc_res*mc_scale_factor, 100 * count_nuebar_cc_res / init_count_nuebar_cc_res);
            printf (" %-12s: %-7.2f(%2.1f%%) %-12s: %-7.2f(%2.1f%%)\n", "Nue CC DIS",   count_nue_cc_dis*mc_scale_factor, 100 * count_nue_cc_dis / init_count_nue_cc_dis, "Nuebar CC DIS",   count_nuebar_cc_dis*mc_scale_factor, 100 * count_nuebar_cc_dis / init_count_nuebar_cc_dis);
            printf (" %-12s: %-7.2f(%2.1f%%) %-12s: %-7.2f(%2.1f%%)\n", "Nue CC COH",   count_nue_cc_coh*mc_scale_factor, 100 * count_nue_cc_coh / init_count_nue_cc_coh, "Nuebar CC COH",   count_nuebar_cc_coh*mc_scale_factor, 100 * count_nuebar_cc_coh / init_count_nuebar_cc_coh);
            printf (" %-12s: %-7.2f(%2.1f%%) %-12s: %-7.2f(%2.1f%%)\n", "Nue CC MEC",   count_nue_cc_mec*mc_scale_factor, 100 * count_nue_cc_mec / init_count_nue_cc_mec, "Nuebar CC MEC",   count_nuebar_cc_mec*mc_scale_factor, 100 * count_nuebar_cc_mec / init_count_nuebar_cc_mec);
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
        
        if (_util.print_mc){
            std::cout << "------------------------------------------------" << std::endl;
            efficiency = double(count_nue_cc + count_nuebar_cc) / double(tot_true_infv_nues);
            purity     = double(count_nue_cc + count_nuebar_cc) / double(sum_mc_dirt_ext);
            printf (" %-20s: ( %-6.1f / %-7.1f ) = %-3.2f %% \n", "Purity", count_nue_cc + count_nuebar_cc, sum_mc_dirt_ext, 100 * purity);
            std::cout << std::endl;
            printf (" %-20s: ( %-6.1f / %-7.1f ) = %-3.2f +/- %-3.2f %% \n", "Efficiency", count_nue_cc + count_nuebar_cc, tot_true_infv_nues, 100 * efficiency, 100 * vec_err.at(k_eff_nu_E_single_bin).at(p));
            printf (" %-20s: ( %-6.1f / %-7.1f ) = %-3.2f +/- %-3.2f %% \n", "Efficiency nue", vec_n.at(k_eff_nu_E_nue_single_bin).at(p), vec_N.at(k_eff_nu_E_nue_single_bin).at(p), 100 * (vec_n.at(k_eff_nu_E_nue_single_bin).at(p) / vec_N.at(k_eff_nu_E_nue_single_bin).at(p)),  100 * vec_err.at(k_eff_nu_E_nue_single_bin).at(p));
            printf (" %-20s: ( %-6.1f / %-7.1f ) = %-3.2f +/- %-3.2f %% \n", "Efficiency nuebar", vec_n.at(k_eff_nu_E_nuebar_single_bin).at(p), vec_N.at(k_eff_nu_E_nuebar_single_bin).at(p), 100 * (vec_n.at(k_eff_nu_E_nuebar_single_bin).at(p) / vec_N.at(k_eff_nu_E_nuebar_single_bin).at(p)), 100 * vec_err.at(k_eff_nu_E_nuebar_single_bin).at(p));
            
            std::cout << std::endl;
            std::cout << " Purity Change     : " <<  100 * (purity - purity_last) << " \%" << std::endl;
            std::cout << " Efficiency Change : "  << 100 * (efficiency - efficiency_last) << " \%" << std::endl;

            efficiency_last = efficiency;
            purity_last     = purity;
        }

        if (_util.print_data){
            std::cout << "------------------------------------------------" << std::endl;
            double tot_bkg = 0;
            if (_util.print_mc && _util.print_dirt && _util.print_ext) tot_bkg = (count_nu_out_fv + count_cosmic + count_numu_cc + count_numu_cc_pi0 + count_nc + count_nc_pi0 + count_unmatched) * mc_scale_factor +
                              count_dirt * dirt_scale_factor + count_ext * ext_scale_factor;
            
            std::cout << " Total Selected Data : " << count_data << std::endl;
            
            if (_util.print_mc && _util.print_dirt && _util.print_ext) std::cout << " Total Nue Candidates in data : " << count_data - tot_bkg<< std::endl;
            
            std::cout << "------------------------------------------------" << std::endl;
        }

        // Fill the efficiency tree
        if (_util.print_mc && _util.print_data && _util.print_ext && _util.print_dirt){
            eff_err = vec_err.at(k_eff_nu_E_single_bin).at(p);
            eff_tree->Fill();
        }

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
        prev_count_thr_nue          = count_thr_nue;
        prev_count_thr_nuebar       = count_thr_nuebar;
        if (_util.print_ext) prev_count_ext = count_ext;
        if (_util.print_dirt) prev_count_dirt = count_dirt;
        
    }

    // Save the efficiency tree
    if (_util.print_mc && _util.print_data && _util.print_ext && _util.print_dirt) {
        
        if (_util.print_mc){
            f_mc->cd();
            eff_tree->Write("",TObject::kOverwrite);
        }
        
    }

}
// -----------------------------------------------------------------------------
void PrintHelper::GetHists(){

    // First get the efficiency histograms
    TEfficiency_hists.resize(k_TH1D_eff_MAX);
    vec_n.resize(k_TH1D_eff_MAX);
    vec_N.resize(k_TH1D_eff_MAX);
    vec_err.resize(k_TH1D_eff_MAX);


    for (unsigned int v = 0; v < TEfficiency_hists.size(); v++){
        TEfficiency_hists.at(v).resize(_util.k_cuts_MAX);
        vec_n.at(v).resize(_util.k_cuts_MAX);
        vec_N.at(v).resize(_util.k_cuts_MAX);
        vec_err.at(v).resize(_util.k_cuts_MAX);
    }
    
    f_mc_hist->cd();
    
    for (unsigned int l = 0; l < _util.k_cuts_MAX; l++ ){
        // Get the efficiency histograms
        _util.GetHist(f_mc_hist, TEfficiency_hists.at(k_eff_nu_E_single_bin).at(l),        Form("TEff/h_true_nu_E_single_bin_%s",_util.cut_dirs.at(l).c_str()));
        _util.GetHist(f_mc_hist, TEfficiency_hists.at(k_eff_nu_E_nue_single_bin).at(l),    Form("TEff/h_true_nu_E_nue_single_bin_%s",_util.cut_dirs.at(l).c_str()));
        _util.GetHist(f_mc_hist, TEfficiency_hists.at(k_eff_nu_E_nuebar_single_bin).at(l), Form("TEff/h_true_nu_E_nuebar_single_bin_%s",_util.cut_dirs.at(l).c_str()));

        // nue + nuebar
        double n, N;
        vec_n.at(k_eff_nu_E_single_bin).at(l)   = TEfficiency_hists.at(k_eff_nu_E_single_bin).at(l)->GetBinContent(1); // total number of signal events for cut i
        vec_N.at(k_eff_nu_E_single_bin).at(l)   = TEfficiency_hists.at(k_eff_nu_E_single_bin).at(0)->GetBinContent(1); // total number of generated events for cut i
        n = vec_n.at(k_eff_nu_E_single_bin).at(l);
        N = vec_N.at(k_eff_nu_E_single_bin).at(l);
        vec_err.at(k_eff_nu_E_single_bin).at(l) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

        // nue
        vec_n.at(k_eff_nu_E_nue_single_bin).at(l)   = TEfficiency_hists.at(k_eff_nu_E_nue_single_bin).at(l)->GetBinContent(1);
        vec_N.at(k_eff_nu_E_nue_single_bin).at(l)   = TEfficiency_hists.at(k_eff_nu_E_nue_single_bin).at(0)->GetBinContent(1);
        n = vec_n.at(k_eff_nu_E_nue_single_bin).at(l);
        N = vec_N.at(k_eff_nu_E_nue_single_bin).at(l);
        vec_err.at(k_eff_nu_E_nue_single_bin).at(l) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

        // nuebar
        vec_n.at(k_eff_nu_E_nuebar_single_bin).at(l)   = TEfficiency_hists.at(k_eff_nu_E_nuebar_single_bin).at(l)->GetBinContent(1);
        vec_N.at(k_eff_nu_E_nuebar_single_bin).at(l)   = TEfficiency_hists.at(k_eff_nu_E_nuebar_single_bin).at(0)->GetBinContent(1);
        n = vec_n.at(k_eff_nu_E_nuebar_single_bin).at(l);
        N = vec_N.at(k_eff_nu_E_nuebar_single_bin).at(l);
        vec_err.at(k_eff_nu_E_nuebar_single_bin).at(l) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

    }

}
// -----------------------------------------------------------------------------
