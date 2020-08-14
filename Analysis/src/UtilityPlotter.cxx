#include "../include/UtilityPlotter.h"

// -----------------------------------------------------------------------------
void UtilityPlotter::Initialise(const char *_run_period, Utility _utility, const char* mode){

    std::cout << "Initalising Utility Plotter ..." << std::endl;
    _util = _utility;

    // Set the run period
    run_period = std::string(_run_period);

    f_nuexsec    = TFile::Open( Form("files/trees/nuexsec_tree_merged_run%s.root", _run_period ));
        
    // Get the Ttree
    _util.GetTree(f_nuexsec, tree, "tree");
    
    // Initialise the TTtree
    InitTree();

    // Standard variation mode
    if (std::string(mode) == "default")  {

        // Look to see if the shower with the most hits is the same as the shower with the most energy
        CompareHitstoEnergy();

        // Lets see how many of the leading showers that we select are not an electorn
        CompareSignalPurity();

        // Make the bin resolution plots
        GetBinResolutions();
        
    }
    // This will call the code to optimise the bin widths
    else if (std::string(mode) == "bins"){
        OptimiseBins();
        return;
    }
    else {
        std::cout << "Error I dont know what mode you have configured..." << mode << std::endl;
        return;
    }
    
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
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareSignalPurity(){

    std::cout <<"\nStudy: Selected shower the electron? \n" << std::endl;

    double total_signal_elec = 0;
    double total_bkg_elec = 0; // Cases where the selected leading shower was not an electron


    // Loop over the entries in the TTree
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 


        // Lets only look at Signal events
        if ((*classifcation == "nue_cc" || *classifcation == "nuebar_cc" || *classifcation == "unmatched_nue" || *classifcation == "unmatched_nuebar") && gen == false) {

            if (shr_bkt_pdg == -11 || shr_bkt_pdg == 11) total_signal_elec++;
            else total_bkg_elec++;

        } // end if signal

    } // End event loop

    std::cout << "\nTotal true elec selected: " << total_signal_elec << std::endl;
    std::cout << "\nTotal true bkg selected: " << total_bkg_elec << std::endl;
    std::cout << "\nPecentage of selected showers that are not an electron: " << 100*total_bkg_elec/(total_bkg_elec+total_signal_elec) << std::endl;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::GetFitResult(double &mean, double &sigma, float bin_lower_edge, float bin_upper_edge, TTree* tree, bool save_hist, bool &converged, bool draw_fit_results){
    
    TCut generic_query = "(classifcation.c_str()==\"nue_cc\" || classifcation.c_str()==\"nuebar_cc\") && !gen"; // This gets selected signal events
    TCut bin_query = Form("shr_energy_cali > %f && shr_energy_cali < %f", bin_lower_edge, bin_upper_edge);
    
    TCanvas * c = new TCanvas(Form("c_%f_%f", bin_upper_edge, sigma), "c", 500, 500);

    // Draw the Query
    tree->Draw("elec_e", generic_query && bin_query);
    
    // Get the histogram from the pad
    TH1D *htemp = (TH1D*)gPad->GetPrimitive("htemp");
    
    // Fit it with a Gaussian
    htemp->Fit("gaus");
    
    // Get the fit result
    TF1 *fit_gaus = htemp->GetFunction("gaus");
    
    // Draw the histogram
    if (save_hist) {
        htemp->SetLineWidth(2);
        htemp->SetLineColor(kBlack);
    }
    htemp->Draw("hist");
    if (converged) fit_gaus->Draw("same");

    mean  = fit_gaus->GetParameter(1);
	sigma = fit_gaus->GetParameter(2);

    if (sigma*2 >= bin_upper_edge - bin_lower_edge - 0.01 && sigma*2 <= bin_upper_edge - bin_lower_edge + 0.01){
        std::cout << "Fit has converged!: " << 2*sigma/(bin_upper_edge - bin_lower_edge) << std::endl;
        converged = true;
    }

    TLatex* range;
    TLatex* fit_params;
    if (save_hist){
        range = new TLatex(0.88,0.86, Form("Reco Energy %0.2f - %0.2f GeV",bin_lower_edge, bin_upper_edge ));
        range->SetTextColor(kGray+2);
        range->SetNDC();
        range->SetTextSize(0.038);
        range->SetTextAlign(32);
        range->Draw();

        fit_params = new TLatex(0.88,0.92, Form("Fit Mean: %0.2f GeV, Fit Sigma: %0.2f GeV",mean, sigma ));
        fit_params->SetTextColor(kGray+2);
        fit_params->SetNDC();
        fit_params->SetTextSize(0.038);
        fit_params->SetTextAlign(32);
        if (draw_fit_results) fit_params->Draw();
        htemp->SetTitle("; Truth Electron Energy [GeV]; Entries");
        htemp->SetStats(kFALSE);
        _util.IncreaseLabelSize(htemp, c);
        c->SetTopMargin(0.11);
        c->Print(Form("plots/Binning/bins_%0.2fGeV_to_%0.2f_GeV.pdf",bin_lower_edge, bin_upper_edge ));
    } 


}
// -----------------------------------------------------------------------------
void UtilityPlotter::OptimiseBins(){

    // Create the Bins directory for saving the plots to
    _util.CreateDirectory("Binning", run_period);

    // Load in the tfile and tree
    double mean{0.0}, sigma{0.0};
    bool converged = false;

    // Were do we want to start the fit iteraction from?
    // Generally choose the first bin width to be 0.25 GeV
    float lower_bin = 0.001;
    // float lower_bin = 1.55;
    
    // Loop over the bins
    for (float bin = 0; bin < 8; bin++ ){
        std::cout << "\n\033[0;34mTrying to optimise the next bin\033[0m\n"<< std::endl;
        converged = false;

        // Slide upper bin value till we get 2xthe STD of the fit
        for (float i = lower_bin+0.1; i <= 4.0; i+=0.001) {
            std::cout << "\n\033[0;34mTrying Bin: " << i << "GeV\033[0m\n"<< std::endl;

            // call function which draws the tree to a canvas, fits the tree and returns the fit parameter
            // If the fit has 2xSTD = the reco bin size then we have successfully optimised the bin
            GetFitResult(mean, sigma, lower_bin, i, tree, false, converged, false);

            // If it converged, do it again and print the canvas then break
            if (converged) {
                GetFitResult(mean, sigma, lower_bin, i, tree, true, converged, true);
                std::cout << "\n\033[0;34mMean: " << mean << "  Sigma: " << sigma<< "\033[0m\n"<< std::endl;
                
                // Reset the lower bin value
                lower_bin = i;
                break;
            }

            // Since the fit doesnt want to converge for the last bin, lets jsut draw it anyway
            if (bin == 7){
                GetFitResult(mean, sigma, 2.63, 3.5, tree, true, converged, false);
                break;
            }

        }
    }
    
}
// -----------------------------------------------------------------------------
void UtilityPlotter::GetBinResolutions(){

    // Create the resolutions directory for saving the plots to
    _util.CreateDirectory("Resolution", run_period);


    // Get the vector of bins
    std::vector<double> bins = _util.reco_shr_bins;
    
    // Loop over the bins
    for (float bin = 0; bin < bins.size()-1; bin++ ){

        std::cout <<"Bin Range: " << bins.at(bin) << " - " << bins.at(bin+1) << " GeV" << std::endl;
        
        PlotResolution(bins.at(bin), bins.at(bin+1), tree, "reco");
        PlotResolution(bins.at(bin), bins.at(bin+1), tree, "true");


    }
    
}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotResolution(float bin_lower_edge, float bin_upper_edge, TTree* tree, std::string variable){
    
    TCut generic_query = "(classifcation.c_str()==\"nue_cc\" || classifcation.c_str()==\"nuebar_cc\") && !gen && elec_e > 0"; // This gets selected events in the MC
    TCut bin_query = Form("shr_energy_cali > %f && shr_energy_cali < %f", bin_lower_edge, bin_upper_edge);
    
    TCanvas * c = new TCanvas(Form("c_%f_%f_%s", bin_upper_edge, bin_lower_edge, variable.c_str()), "c", 500, 500);

    TH1D *htemp = new TH1D("htemp","", 30, -3, 3);

    // Draw the Query
    if (variable == "reco") tree->Draw("shr_energy_cali - elec_e / shr_energy_cali >> htemp", generic_query && bin_query);
    else                    tree->Draw("shr_energy_cali - elec_e / elec_e >> htemp", generic_query && bin_query);
    
    // Get the histogram from the pad
   
            
    // Draw the histogram
    htemp->SetLineWidth(2);
    htemp->SetLineColor(kBlack);
    htemp->Draw("hist");

    TLatex* range;
    range = new TLatex(0.88,0.86, Form("Reco Energy %0.2f - %0.2f GeV",bin_lower_edge, bin_upper_edge ));
    range->SetTextColor(kGray+2);
    range->SetNDC();
    range->SetTextSize(0.038);
    range->SetTextAlign(32);
    range->Draw();

    if (variable == "reco") htemp->SetTitle("; Reco - True / Reco; Entries");
    else                    htemp->SetTitle("; Reco - True / True; Entries");
    htemp->SetStats(kFALSE);
    _util.IncreaseLabelSize(htemp, c);
    c->SetTopMargin(0.11);
    if (variable== "reco") c->Print(Form("plots/run%s/Resolution/resolution_%0.2fGeV_to_%0.2f_GeV_reco.pdf", run_period.c_str(), bin_lower_edge, bin_upper_edge ));
    else c->Print(Form("plots/run%s/Resolution/resolution_%0.2fGeV_to_%0.2f_GeV_true.pdf", run_period.c_str(), bin_lower_edge, bin_upper_edge ));
    
    delete htemp;


}
// -----------------------------------------------------------------------------