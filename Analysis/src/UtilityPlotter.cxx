#include "../include/UtilityPlotter.h"

// -----------------------------------------------------------------------------
void UtilityPlotter::Initialise(const char *_run_period, Utility _utility, const char* mode){

    std::cout << "Initalising Utility Plotter ..." << std::endl;
    _util = _utility;

    // Set the run period
    run_period = std::string(_run_period);

    // Standard variation mode
    if (std::string(mode) == "default")  {
        f_nuexsec    = TFile::Open( Form("files/trees/nuexsec_tree_merged_run%s.root", _run_period ));
        
        // Get the Ttree
        _util.GetTree(f_nuexsec, tree, "tree");
        
        // Initialise the TTtree
        InitTree();
    }
    else {
        std::cout << "Error I dont know what mode you have configured..." << mode << std::endl;
        exit(1);
    }
    
    // Now we have initialised the Ttree, we can do some studies!
    
    // Look to see if the shower with the most hits is the same as the shower with the most energy
    CompareHitstoEnergy();

}
// -----------------------------------------------------------------------------
void UtilityPlotter::InitTree(){
    
    // Set the tree branches
    tree->SetBranchAddress("gen",    &gen);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("true_energy", &true_energy);
    tree->SetBranchAddress("reco_energy", &reco_energy);
    tree->SetBranchAddress("classifcation",   &classifcation);
    tree->SetBranchAddress("shr_energy_cali", &shr_energy_cali);
    tree->SetBranchAddress("elec_e",  &elec_e);
    tree->SetBranchAddress("ppfx_cv",  &ppfx_cv);
    tree->SetBranchAddress("weightSplineTimesTune",  &weightSplineTimesTune);
    tree->SetBranchAddress("numi_ang",  &numi_ang);
    tree->SetBranchAddress("nu_pdg",  &nu_pdg);
    tree->SetBranchAddress("shr_bkt_purity",         &shr_bkt_purity);
    tree->SetBranchAddress("shr_bkt_completeness",   &shr_bkt_completeness);
    tree->SetBranchAddress("shr_bkt_E",              &shr_bkt_E);
    tree->SetBranchAddress("shr_bkt_pdg", &shr_bkt_pdg);
    tree->SetBranchAddress("all_shr_hits" , &all_shr_hits);
    tree->SetBranchAddress("all_shr_energies",  &all_shr_energies);
}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareHitstoEnergy(){

    std::cout <<"\nStudy: Tot hits to leading shower? \n" << std::endl;

    double total_signal_events = 0;
    double total_non_leading_hit_events_sig = 0; // The total number of events where the most number of hits is the most energetic shower

    double total_bkg_events = 0;
    double total_non_leading_hit_events_bkg = 0;

    double total_events = 0;
    double total_non_leading_hit_events = 0;

    // Loop over the entries in the TTree
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 


        // Lets only look at Signal events
        if ((*classifcation == "nue_cc" || *classifcation == "nuebar_cc" || *classifcation == "unmatched_nue" || *classifcation == "unmatched_nuebar") && gen == false) {

            // Get the leading shower index
            double index_leading_shr_hits = 0;
            double leading_shr_hits = 0;
            double leading_energy_hits = 0;
            
            double index_leading_shr_E = 0;
            double leading_shr_E = 0;

            if (all_shr_hits->size() != all_shr_energies->size()) std::cout <<"Warning hit vector size does not equal shower energy vector size!" <<std::endl;

            for (unsigned int index = 0 ; index < all_shr_hits->size(); index++){

                // Set the leading shower
                if (all_shr_hits->at(index) > leading_shr_hits){
                    leading_shr_hits = all_shr_hits->at(index);
                    index_leading_shr_hits = index;
                }

                // Set the leading shower energy
                if (all_shr_energies->at(index) > leading_shr_E){
                    leading_shr_E = all_shr_energies->at(index);
                    leading_energy_hits = all_shr_hits->at(index);
                    index_leading_shr_E = index;
                }

            }

            if (index_leading_shr_E != index_leading_shr_hits){
                std::cout << *classifcation << " Shr Hits: " << leading_shr_hits << "  Shr Energy: " << leading_shr_E <<  "  Leading energy hits: " << leading_energy_hits<< std::endl;
                total_non_leading_hit_events_sig++;
                total_non_leading_hit_events++;
            }


            total_signal_events++;
            total_events++;

        } // end if signal

        // Background event
        if ( *classifcation == "nu_out_fv"  || *classifcation == "cosmic"      ||
                *classifcation == "numu_cc"    || *classifcation == "numu_cc_pi0" || *classifcation == "nc" || 
                *classifcation == "nc_pi0"     || *classifcation == "cosmic_nue" || *classifcation == "cosmic_nuebar"){
            
            // Get the leading shower index
            double index_leading_shr_hits = 0;
            double leading_shr_hits = 0;
            double leading_energy_hits = 0;
            
            double index_leading_shr_E = 0;
            double leading_shr_E = 0;

            if (all_shr_hits->size() != all_shr_energies->size()) std::cout <<"Warning hit vector size does not equal shower energy vector size!" <<std::endl;

            for (unsigned int index = 0 ; index < all_shr_hits->size(); index++){

                // Set the leading shower
                if (all_shr_hits->at(index) > leading_shr_hits){
                    leading_shr_hits = all_shr_hits->at(index);
                    index_leading_shr_hits = index;
                }

                // Set the leading shower energy
                if (all_shr_energies->at(index) > leading_shr_E){
                    leading_shr_E = all_shr_energies->at(index);
                    leading_energy_hits = all_shr_hits->at(index);
                    index_leading_shr_E = index;
                }

            }

            if (index_leading_shr_E != index_leading_shr_hits){
                std::cout << *classifcation <<  " Shr Hits: " << leading_shr_hits << "  Shr Energy: " << leading_shr_E <<  "  Leading energy hits: " << leading_energy_hits<< std::endl;
                total_non_leading_hit_events_bkg++;
                total_non_leading_hit_events++;
            }


            total_bkg_events++;
            total_events++;
            
        }// End if background event

        

    } // End event loop

    std::cout << "\nPercentage of signal events where shower with most hits is not the most energetic shower: " << 100*total_non_leading_hit_events_sig / total_signal_events << std::endl;
    std::cout << "\nPercentage of background events where shower with most hits is not the most energetic shower: " << 100*total_non_leading_hit_events_bkg / total_bkg_events << std::endl;
    std::cout << "\nPercentage of all events where shower with most hits is not the most energetic shower: " << 100*total_non_leading_hit_events / total_events << std::endl;

}
