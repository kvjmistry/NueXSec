#include "../include/PrintHelper.h"

// -----------------------------------------------------------------------------
void PrintHelper::Initialise(const char* run_period, const char * mc_file_in, bool _print_mc, bool _print_data, bool _print_ext, bool _print_dirt, utility _utility ){

    std::cout << "Initalising Print Helper..." << std::endl;

    _util = _utility;

    std::string file_out_str = mc_file_in;

    std::string file_name;

    print_mc    = _print_mc;
    print_data  = _print_data;
    print_ext   = _print_ext;
    print_dirt  = _print_dirt;

    // Set the scale factors
    if (strcmp(run_period, "1") == 0){
        mc_scale_factor     = _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_MC_POT);
        dirt_scale_factor   = _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_Dirt_POT);
        intime_scale_factor = _util.config_v.at(_util.k_Run1_Data_trig) / _util.config_v.at(_util.k_Run1_EXT_trig);
    }
    else if (strcmp(run_period, "3") == 0){
        mc_scale_factor     = _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_MC_POT);
        dirt_scale_factor   = _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_Dirt_POT);
        intime_scale_factor = _util.config_v.at(_util.k_Run3_Data_trig) / _util.config_v.at(_util.k_Run3_EXT_trig);
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
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
        
        mc_counter_tree->SetBranchAddress("count_nue_cc",       &count_nue_cc);
        mc_counter_tree->SetBranchAddress("count_nue_cc_mixed", &count_nue_cc_mixed);
        mc_counter_tree->SetBranchAddress("count_nu_out_fv",    &count_nu_out_fv);
        mc_counter_tree->SetBranchAddress("count_cosmic",       &count_cosmic);
        mc_counter_tree->SetBranchAddress("count_numu_cc",      &count_numu_cc);
        mc_counter_tree->SetBranchAddress("count_numu_cc_pi0",  &count_numu_cc_pi0);
        mc_counter_tree->SetBranchAddress("count_nc",           &count_nc);
        mc_counter_tree->SetBranchAddress("count_nc_pi0",       &count_nc_pi0);
        mc_counter_tree->SetBranchAddress("count_unmatched",    &count_unmatched);
        mc_counter_tree->SetBranchAddress("count_total_mc",     &count_total_mc);
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

        if (p == 0) tot_true_infv_nues = count_nue_cc + count_nue_cc_mixed;

        
        // // Sum of selected mc, dirt and ext. The dirt and ext are scaled to the MC POT
        double sum_mc_dirt_ext = count_total_mc+ (count_ext * (intime_scale_factor / mc_scale_factor)) + (count_dirt * (dirt_scale_factor / mc_scale_factor));

        if (print_mc && print_data && print_ext && print_dirt){
            std::cout << "\n------------------------------------------------" << std::endl;
        }
        std::cout << "------------------------------------------------" << std::endl;
        
        std::cout.precision(6);
        std::cout << "\n\033[0;33m <" << _util.cut_dirs.at(p) << "> \033[0m" << std::endl;
        
        if (!print_mc) {
            if (p == 0) std::cout << "                    Scaled to MC POT | Scaled to Data POT | Unscaled" << std::endl;
        }
        else std::cout << "                    Scaled to MC POT | Scaled to Data POT | Unscaled" << std::endl;

        if (print_mc){
            std::cout << " Total Candidate Nue     : " << sum_mc_dirt_ext                 << "    " << double(sum_mc_dirt_ext                * mc_scale_factor  ) << std::endl;
            std::cout << " Number of Nue CC        : " << count_nue_cc          << "    " << double(count_nue_cc         * mc_scale_factor  ) << std::endl;
            std::cout << " Number of Nue CC Mixed  : " << count_nue_cc_mixed    << "    " << double(count_nue_cc_mixed   * mc_scale_factor  ) << std::endl;
            std::cout << " Number of Nu out FV     : " << count_nu_out_fv       << "    " << double(count_nu_out_fv      * mc_scale_factor  ) << std::endl;
            std::cout << " Number of Cosmic        : " << count_cosmic          << "    " << double(count_cosmic         * mc_scale_factor  ) << std::endl;
            std::cout << " Number of Numu CC       : " << count_numu_cc         << "    " << double(count_numu_cc        * mc_scale_factor  ) << std::endl;
            std::cout << " Number of Numu CC Pi0   : " << count_numu_cc_pi0     << "    " << double(count_numu_cc_pi0    * mc_scale_factor  ) << std::endl;
            std::cout << " Number of NC            : " << count_nc              << "    " << double(count_nc             * mc_scale_factor  ) << std::endl;
            std::cout << " Number of NC Pi0        : " << count_nc_pi0          << "    " << double(count_nc_pi0         * mc_scale_factor  ) << std::endl;
            std::cout << " Number of Unmatched     : " << count_unmatched       << "    " << double(count_unmatched      * mc_scale_factor  ) << std::endl;
        }

        if (print_ext) { 
            std::cout << " Number of InTime Cosmics: " << double(count_ext * (intime_scale_factor / mc_scale_factor))
                << "\t " << double(count_ext * intime_scale_factor) << "\t " << count_ext << std::endl;
        }
        
        if (print_dirt){
            std::cout << " Number of Dirt          : " << double(count_dirt * dirt_scale_factor / mc_scale_factor)
                << "\t "                         << double(count_dirt * dirt_scale_factor) << "\t " << count_dirt << std::endl;
        }
        
        if (print_mc){
            std::cout << "----------- Neutrinos in FV Truth -------------" << std::endl;
            std::cout << " Nue CC QE    : " << count_nue_cc_qe   <<  "   Nuebar CC QE    : " << count_nuebar_cc_qe  << std::endl;
            std::cout << " Nue CC Res   : " << count_nue_cc_res  <<  "   Nuebar CC Res   : " << count_nuebar_cc_res << std::endl;
            std::cout << " Nue CC DIS   : " << count_nue_cc_dis  <<  "   Nuebar CC DIS   : " << count_nuebar_cc_dis << std::endl;
            std::cout << " Nue CC COH   : " << count_nue_cc_coh  <<  "   Nuebar CC COH   : " << count_nuebar_cc_coh << std::endl;
            std::cout << " Nue CC MEC   : " << count_nue_cc_mec  <<  "   Nuebar CC MEC   : " << count_nuebar_cc_mec << std::endl;
            std::cout << std::endl;
            std::cout << " Numu CC QE   : " << count_numu_cc_qe  <<  "   Numubar CC QE   : " << count_numubar_cc_qe << std::endl;
            std::cout << " Numu CC Res  : " << count_numu_cc_res <<  "   Numubar CC Res  : " << count_numubar_cc_res<< std::endl;
            std::cout << " Numu CC DIS  : " << count_numu_cc_dis <<  "   Numubar CC DIS  : " << count_numubar_cc_dis<< std::endl;
            std::cout << " Numu CC COH  : " << count_numu_cc_coh <<  "   Numubar CC COH  : " << count_numubar_cc_coh<< std::endl;
            std::cout << " Numu CC MEC  : " << count_numu_cc_mec <<  "   Numubar CC MEC  : " << count_numubar_cc_mec<< std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << " Tot Nue in FV                 : " << count_nue_cc_infv     << "   Tot Nue in Cryo                 : " << count_nue_cc_incryo     << std::endl;
            std::cout << " Tot Nuebar in FV              : " << count_nuebar_cc_infv  << "   Tot Nuebar in Cryo              : " << count_nuebar_cc_incryo << std::endl;
            std::cout << " Tot NuMu in FV                : " << count_numu_cc_infv    << "   Tot NuMu in Cryo                : " << count_numu_cc_incryo<< std::endl;
            std::cout << " Tot NuMubar in FV             : " << count_numubar_cc_infv << "   Tot NuMubar in Cryo             : " << count_numubar_cc_incryo<< std::endl;
            // std::cout << " Tot NC                  : " << counter_tot_nue_numu_nc << std::endl;
            // std::cout << " Sum Neutrinos           : " <<  << std::endl;
        }
        
        if (print_mc){
            std::cout << "------------------------------------------------" << std::endl;
            efficiency = double(count_nue_cc + count_nue_cc_mixed) / double(tot_true_infv_nues);
            purity     = double(count_nue_cc + count_nue_cc_mixed) / double(sum_mc_dirt_ext);
            std::cout << " Efficiency       : " << "( " << count_nue_cc + count_nue_cc_mixed << " / " << tot_true_infv_nues << " ) = " << efficiency << std::endl;
            std::cout << " Purity           : " << "( " << count_nue_cc + count_nue_cc_mixed << " / " << sum_mc_dirt_ext           << " ) = " << purity << std::endl;
        }

        if (print_data){
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << " Total Nue Candidates in data : " << count_data << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
        }

        // Fill the efficiency tree
        if (print_mc && print_data && print_ext && print_dirt) eff_tree->Fill();
        
    }

    // Save the efficiency tree
    if (print_mc && print_data && print_ext && print_dirt) {
        
        if (print_mc){
            f_mc->cd();
            eff_tree->Write("",TObject::kOverwrite);
        }
        
    }

}
