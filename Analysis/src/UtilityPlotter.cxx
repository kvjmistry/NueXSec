#include "../include/UtilityPlotter.h"

// -----------------------------------------------------------------------------
void UtilityPlotter::Initialise(Utility _utility){

    std::cout << "Initalising Utility Plotter ..." << std::endl;
    _util = _utility;

    f_nuexsec    = TFile::Open( Form("files/trees/nuexsec_tree_merged_run%s.root", _util.run_period ));
        
    // Get the Ttree
    _util.GetTree(f_nuexsec, tree, "tree");
    
    // Initialise the TTtree
    InitTree();

    // Standard variation mode
    if (std::string(_util.uplotmode) == "default")  {

        // Compare the efficiency for the det var CV and intrinsic nue det var CV
        // leave this commented out, needs speccific files for this to run 
        // CompareDetVarEfficiency();

        // Compare the efficiency in run 1 and run 3
        CompareEfficiency();

        // Function that plots all the ppfx universe weights on one plot for the backgrounds
        StudyPPFXWeights();

        // Look to see if the shower with the most hits is the same as the shower with the most energy
        CompareHitstoEnergy();

        // Lets see how many of the leading showers that we select are not an electorn
        CompareSignalPurity();

        // Make the bin resolution plots
        PlotVarbyRecoBin();

        // Plot the 1D flux with the threhsold line
        PlotIntegratedFluxwithThrehold();
        
    }
    // Make the true variable plots
    else if (std::string(_util.uplotmode) == "true"){
        PlotTrueVar();
        return;

    }
    // This will call the code to optimise the bin widths
    else if (std::string(_util.uplotmode) == "bins"){
        OptimiseBins();
        return;
    }
    // This will call the code to optimise the bin widths
    else if (std::string(_util.uplotmode) == "models"){ 
        // Set the names of the histograms
        _util.SetAxesNames(var_labels_xsec, var_labels_events, var_labels_eff, smear_hist_name, vars, xsec_scale);

        _util.CreateDirectory("Models/" + std::string(_util.xsec_var));
        TestModelDependence();
        CompareDataCrossSections();
        CompareSmearing();
        CompareUnfoldedModels();
        return;
    }
    else {
        std::cout << "Error I dont know what mode you have configured..." << _util.uplotmode << std::endl;
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
    tree->SetBranchAddress("classification",   &classification);
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
    tree->SetBranchAddress("weightsPPFX",           &weightsPPFX);
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
        if ((*classification == "nue_cc" || *classification == "nuebar_cc" || *classification == "unmatched_nue" || *classification == "unmatched_nuebar") && gen == false) {

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
                std::cout << *classification << " Shr Hits: " << leading_shr_hits << "  Shr Energy: " << leading_shr_E <<  "  Leading energy hits: " << leading_energy_hits<< std::endl;
                total_non_leading_hit_events_sig++;
                total_non_leading_hit_events++;
            }


            total_signal_events++;
            total_events++;

        } // end if signal

        // Background event
        if ( *classification == "nu_out_fv"  || *classification == "cosmic"      ||
                *classification == "numu_cc"    || *classification == "numu_cc_pi0" || *classification == "nc" || 
                *classification == "nc_pi0"     || *classification == "cosmic_nue" || *classification == "cosmic_nuebar"){
            
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
                std::cout << *classification <<  " Shr Hits: " << leading_shr_hits << "  Shr Energy: " << leading_shr_E <<  "  Leading energy hits: " << leading_energy_hits<< std::endl;
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
        if ((*classification == "nue_cc" || *classification == "nuebar_cc" || *classification == "unmatched_nue" || *classification == "unmatched_nuebar") && gen == false) {

            if (shr_bkt_pdg == -11 || shr_bkt_pdg == 11) total_signal_elec++;
            else total_bkg_elec++;

        } // end if signal

    } // End event loop

    std::cout << "\nTotal true elec selected: " << total_signal_elec << std::endl;
    std::cout << "\nTotal true bkg selected: " << total_bkg_elec << std::endl;
    std::cout << "\nPecentage of selected showers that are not an electron: " << 100*total_bkg_elec/(total_bkg_elec+total_signal_elec) << std::endl;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::GetFitResult(double &mean, double &sigma, float bin_lower_edge, float bin_upper_edge, TTree* tree, bool save_hist, bool &converged, bool draw_fit_results, std::string var){
        
    TCut query;
    
    if (var == "elec_E"){
        query = Form(" weight*(  ((classification.c_str()==\"nue_cc\" || classification.c_str()==\"nuebar_cc\") && passed_selection) && shr_energy_cali > %f && shr_energy_cali < %f )", bin_lower_edge, bin_upper_edge); 
    }
    else if (var == "elec_ang"){
        query = Form(" weight*(  ((classification.c_str()==\"nue_cc\" || classification.c_str()==\"nuebar_cc\") && passed_selection) && effective_angle > %f && effective_angle < %f )", bin_lower_edge, bin_upper_edge); 
    }
    else if (var == "elec_cang"){
        query = Form(" weight*(  ((classification.c_str()==\"nue_cc\" || classification.c_str()==\"nuebar_cc\") && passed_selection) && cos_effective_angle > %f && cos_effective_angle < %f )", bin_lower_edge, bin_upper_edge); 
    }
    else {
        return;
    }


    TCanvas * c = new TCanvas(Form("c_%f_%f", bin_upper_edge, sigma), "c", 500, 500);

    // Get the histogram from the pad
    TH1D *htemp;
    
    if (var == "elec_E"){
        htemp = new TH1D("htemp","", 90, 0, 6.0); // Set the binnning
    }
    else if (var == "elec_ang"){
        htemp = new TH1D("htemp","", 90, 0, 180); // Set the binnning
    }
    else if (var == "elec_cang"){
        htemp = new TH1D("htemp","", 90, -1.0, 1.0); // Set the binnning
    }
    else {
        return;
    }

    // Draw the Query and put into histogram
    if (var == "elec_E"){
        tree->Draw("elec_e >> htemp", query);
    }
    else if (var == "elec_ang"){
        tree->Draw("true_effective_angle >> htemp", query);
    }
    else if (var == "elec_cang"){
        tree->Draw("cos_true_effective_angle >> htemp", query);
    }

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


    if (var == "elec_E"){
        if ( (sigma*2 >= bin_upper_edge - bin_lower_edge - 0.01 && sigma*2 <= bin_upper_edge - bin_lower_edge + 0.01) && htemp->Integral() > 330){
            std::cout << "Fit has converged!: " << 2*sigma/(bin_upper_edge - bin_lower_edge) << std::endl;
            converged = true;
        }
    }
    else if (var == "elec_ang"){
        if ( (sigma*2 >= bin_upper_edge - bin_lower_edge - 5.0 && sigma*2 <= bin_upper_edge - bin_lower_edge + 5.0) && htemp->Integral() > 330){
            std::cout << "Fit has converged!: " << 2*sigma/(bin_upper_edge - bin_lower_edge) << std::endl;
            converged = true;
        }
    }
    else if (var == "elec_cang"){
        if ( (sigma*2 >= bin_upper_edge - bin_lower_edge - 0.6 && sigma*2 <= bin_upper_edge - bin_lower_edge + 0.6) && htemp->Integral() > 330){
            std::cout << "Fit has converged!: " << 2*sigma/(bin_upper_edge - bin_lower_edge) << std::endl;
            converged = true;
        }
    }

    TLatex* range;
    TLatex* fit_params;
    if (save_hist){
        
        if (var == "elec_E"){
            range = new TLatex(0.88,0.86, Form("Reco Energy %0.2f - %0.2f GeV",bin_lower_edge, bin_upper_edge ));
        }
        else if (var == "elec_ang"){
            range = new TLatex(0.88,0.86, Form("Reco #beta %0.2f - %0.2f deg",bin_lower_edge, bin_upper_edge ));
        }
         else if (var == "elec_cang"){
            range = new TLatex(0.88,0.86, Form("Reco cos(#beta) %0.2f - %0.2f",bin_lower_edge, bin_upper_edge ));
        }
        
        range->SetTextColor(kGray+2);
        range->SetNDC();
        range->SetTextSize(0.038);
        range->SetTextAlign(32);
        range->Draw();

        if (var == "elec_E"){
            fit_params = new TLatex(0.88,0.92, Form("Fit Mean: %0.2f GeV, Fit Sigma: %0.2f GeV",mean, sigma ));
        }
        else if (var == "elec_ang"){
            fit_params = new TLatex(0.88,0.92, Form("Fit Mean: %0.2f deg, Fit Sigma: %0.2f deg",mean, sigma ));
        }
        else if (var == "elec_cang"){
            fit_params = new TLatex(0.88,0.92, Form("Fit Mean: %0.2f, Fit Sigma: %0.2f",mean, sigma ));
        }

        fit_params->SetTextColor(kGray+2);
        fit_params->SetNDC();
        fit_params->SetTextSize(0.038);
        fit_params->SetTextAlign(32);
        if (draw_fit_results) fit_params->Draw();
        

        if (var == "elec_E"){
            htemp->SetTitle("; E^{true}_{e#lower[-0.5]{-} + e^{+}} [GeV]; Entries");
        }
        else if (var == "elec_ang"){
            htemp->SetTitle("; #beta^{true}_{e#lower[-0.5]{-} + e^{+}} [deg]; Entries");
        }
        else if (var == "elec_cang"){
            htemp->SetTitle("; cos(#beta)^{true}_{e#lower[-0.5]{-} + e^{+}}; Entries");
        }
        
        
        htemp->SetStats(kFALSE);
        _util.IncreaseLabelSize(htemp, c);
        c->SetTopMargin(0.11);
        
        if (var == "elec_E"){
            c->Print(Form("plots/run%s/Binning/%s/bins_%0.2fGeV_to_%0.2fGeV.pdf",_util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        else if (var == "elec_ang"){
            c->Print(Form("plots/run%s/Binning/%s/bins_%0.1fdeg_to_%0.1fdeg.pdf",_util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        else if (var == "elec_cang"){
            c->Print(Form("plots/run%s/Binning/%s/bins_%0.2f_to_%0.2f.pdf",_util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
    } 

    std::cout << _util.red << "Integral: " <<htemp->Integral()  << _util.reset<< std::endl;

    delete htemp;
    delete c;


}
// -----------------------------------------------------------------------------
void UtilityPlotter::OptimiseBins(){

    // Create the Bins directory for saving the plots to
    _util.CreateDirectory("Binning/" + std::string(_util.xsec_var));

    // The variable to optimise
    std::string var = std::string(_util.xsec_var);
    int nbins;

    // Load in the tfile and tree
    double mean{0.0}, sigma{0.0};
    bool converged = false;

    // Were do we want to start the fit iteraction from?
    // Generally choose the first bin width to be 0.25 GeV
    float lower_bin = 0.001;
    // float lower_bin = 0.3;

    // The last bin to go up to
    float upper_bin = 6.0;

    // What increment size to increase the bins by
    float increment_size = 0.001;

    if (var == "elec_E"){
        nbins = 6;
        lower_bin = 0.001;
        upper_bin = 6.0;
        increment_size = 0.001;
    }
    else if (var == "elec_ang"){
        nbins = 6;
        lower_bin = 0.0;
        upper_bin = 180.;
        increment_size = 0.5;
    }
     else if (var == "elec_cang"){
        nbins = 5;
        lower_bin = -1.0;
        upper_bin = 1.0;
        increment_size = 0.005;
    }
    
    // Loop over the bins
    for (float bin = 0; bin < nbins; bin++ ){
        std::cout << "\n\033[0;34mTrying to optimise the next bin\033[0m\n"<< std::endl;
        converged = false;

        // Slide upper bin value till we get 2xthe STD of the fit
        for (float i = lower_bin+increment_size; i <= upper_bin; i+=increment_size) {
            
            if (var == "elec_E")
                std::cout << "\n\033[0;34mTrying Bin: " << i << "GeV\033[0m\n"<< std::endl;
            else if (var == "elec_ang")
                std::cout << "\n\033[0;34mTrying Bin: " << i << "deg\033[0m\n"<< std::endl;
            else if (var == "elec_cang")
                std::cout << "\n\033[0;34mTrying Bin: " << i << "\033[0m\n"<< std::endl;
            else
                std::cout << "Warning, unknown variable specified" << std::endl;

            // Since in the case of electron energy, the first bin does not want to converge
            // We set this manually
            if (bin == 0 && var == "elec_E"){
                bool fake = true;
                GetFitResult(mean, sigma, 0.0, 0.30, tree, true, fake, true, var);
                lower_bin = 0.30;
                break;
            }

            if (bin == 0 && var == "elec_ang"){
                bool fake = true;
                GetFitResult(mean, sigma, 0.0, 6.0, tree, true, fake, true, var);
                lower_bin = 6.0;
                break;
            }

            if (bin == 0 && var == "elec_cang"){
                // bool fake = true;
                GetFitResult(mean, sigma, -1.0, 0.6, tree, true, converged, true, var);
                lower_bin = 0.6;
                break;
            }
            
            // call function which draws the tree to a canvas, fits the tree and returns the fit parameter
            // If the fit has 2xSTD = the reco bin size then we have successfully optimised the bin
            GetFitResult(mean, sigma, lower_bin, i, tree, false, converged, false, var);

            // If it converged, do it again and print the canvas then break
            if (converged) {
                
                GetFitResult(mean, sigma, lower_bin, i, tree, true, converged, true, var);
                std::cout << "\n\033[0;34mMean: " << mean << "  Sigma: " << sigma<< "\033[0m\n"<< std::endl;
                
                // Reset the lower bin value
                lower_bin = i;
                break;
            }

            // Since the fit doesnt want to converge for elec_E for the last bin, lets just draw it anyway
            if (bin == nbins-1){
                bool fake = true;
                GetFitResult(mean, sigma, lower_bin, upper_bin, tree, true, fake, true, var);
                break;
            }

        }
    }
    
}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotVarbyRecoBin(){

    // Create the resolutions directory for saving the plots to
    _util.CreateDirectory("Resolution");

    _util.CreateDirectory("Purity_Completeness");


    // Get the vector of bins
    std::vector<double> bins = _util.reco_shr_bins;
    
    // Loop over the bins
    for (float bin = 0; bin < bins.size()-1; bin++ ){

        std::cout <<"\nBin Range: " << bins.at(bin) << " - " << bins.at(bin+1) << " GeV" << std::endl;
        
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "reco_e");
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "true_e");
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "purity");
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "completeness");


    }
    
}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotQuery(float bin_lower_edge, float bin_upper_edge, TTree* tree, std::string variable_str){
    
    TCut generic_query = "(classification.c_str()==\"nue_cc\" || classification.c_str()==\"nuebar_cc\") && !gen && elec_e > 0"; // This gets selected signal events in the MC
    TCut bin_query = Form("shr_energy_cali > %f && shr_energy_cali < %f", bin_lower_edge, bin_upper_edge); // Get the reconstructed shower energy range

    TCut query = Form(" weight*(  ((classification.c_str()==\"nue_cc\" || classification.c_str()==\"nuebar_cc\") && passed_selection) && shr_energy_cali > %f && shr_energy_cali < %f )", bin_lower_edge, bin_upper_edge); 
    
    TCanvas * c = new TCanvas(Form("c_%f_%f_%s", bin_upper_edge, bin_lower_edge, variable_str.c_str()), "c", 500, 500);

    TH1D *htemp;
    if      (variable_str == "reco_e") htemp = new TH1D("htemp","", 30, -1.2, 1.2);
    else if (variable_str == "true_e") htemp = new TH1D("htemp","", 30, -1.2, 1.2);
    else if (variable_str == "purity") htemp = new TH1D("htemp","", 21, 0, 1.1);
    else if (variable_str == "completeness") htemp = new TH1D("htemp","", 21, 0, 1.1);
    else {
        std::cout << "incorrect variable input" << std::endl;
        return;
    }
     

    // Draw the Query -- adjust by query type
    if      (variable_str == "reco_e") tree->Draw("(shr_energy_cali - elec_e) / shr_energy_cali >> htemp", query);
    else if (variable_str == "true_e") tree->Draw("(shr_energy_cali - elec_e) / elec_e >> htemp", query);
    else if (variable_str == "purity") tree->Draw("shr_bkt_purity >> htemp", query);
    else if (variable_str == "completeness") tree->Draw("shr_bkt_completeness >> htemp", query);
    else {
        std::cout << "incorrect variable input" << std::endl;
        return;
    }
            
    // Draw the histogram
    htemp->SetLineWidth(2);
    htemp->SetLineColor(kBlack);
    htemp->Draw("hist");

    // Draw the text specifying the bin range
    TLatex* range = new TLatex(0.65,0.91, Form("Reco Energy %0.2f - %0.2f GeV",bin_lower_edge, bin_upper_edge ));
    _util.SetTextProperties(range);
    range->Draw();

    double entries = htemp->Integral();
    TLatex* text_entries = new TLatex(0.39,0.86, Form("Entries %4.0f", entries ));
    _util.SetTextProperties(text_entries);
    text_entries->Draw();

    Double_t mean = htemp->GetMean();
    TLatex* text_mean = new TLatex(0.39,0.82, Form("Mean %4.2f", mean));
    _util.SetTextProperties(text_mean);
    text_mean->Draw();

    Double_t rms= htemp->GetRMS();
    TLatex* text_rms = new TLatex(0.39,0.78, Form("STD %4.2f", rms));
    _util.SetTextProperties(text_rms);
    text_rms->Draw();

    if (variable_str == "reco_e")      htemp->SetTitle("; Reco - True / Reco; Entries");
    else if (variable_str == "true_e") htemp->SetTitle("; Reco - True / True; Entries");
    else if (variable_str == "purity") htemp->SetTitle("; Reco Shower Purity; Entries");
    else if (variable_str == "completeness") htemp->SetTitle("; Reco Shower Completeness; Entries");
    else {
        std::cout << "incorrect variable input" << std::endl;
        return;
    }
    htemp->SetStats(kFALSE);
    _util.IncreaseLabelSize(htemp, c);
    c->SetTopMargin(0.11);

    
    
    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);
    

    // Save it 
    if (variable_str== "reco_e")       c->Print(Form("plots/run%s/Resolution/resolution_%0.0fMeV_to_%0.0f_MeV_reco.pdf", _util.run_period, bin_lower_edge*1000, bin_upper_edge*1000 ));
    else if (variable_str == "true_e") c->Print(Form("plots/run%s/Resolution/resolution_%0.0fMeV_to_%0.0f_MeV_true.pdf", _util.run_period, bin_lower_edge*1000, bin_upper_edge*1000 ));
    else if (variable_str == "purity") c->Print(Form("plots/run%s/Purity_Completeness/purity_%0.0fMeV_to_%0.0f_MeV.pdf", _util.run_period, bin_lower_edge*1000, bin_upper_edge*1000 ));
    else if (variable_str == "completeness") c->Print(Form("plots/run%s/Purity_Completeness/completeness_%0.0fMeV_to_%0.0f_MeV.pdf", _util.run_period, bin_lower_edge*1000, bin_upper_edge*1000 ));
    else {
        std::cout << "incorrect variable input" << std::endl;
        return;
    }

    delete htemp;


}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotIntegratedFluxwithThrehold(){
    
    // This changes the plot to average flux or not
    bool draw_averge = true;

    gStyle->SetOptStat(0);


    TFile *f;

    TH1D *h_nue, *h_nuebar;

    double flux_scale_factor{1.0e-4}; // unit conversion of flux from m2 to cm2

    // Get the flux File
    // Switch the file path depending on whether we are on the gpvm or not
    // Initialise the Flux file
    
    std::string flux_file_name;
    if (std::string(_util.run_period) == "1"){
        
        // Switch the file path depending on whether we are on the gpvm or not
        if (!_util.use_gpvm)
            flux_file_name = "Systematics/output_fhc_uboone_run0.root";
        else
            flux_file_name = "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_fhc_run0_merged.root";


        std::cout << "Using Flux file name: \033[0;31m" << flux_file_name << "\033[0m" <<  std::endl;
        bool boolfile = _util.GetFile(f, flux_file_name);
        if (boolfile == false) gSystem->Exit(0); 
    }
    else if (std::string(_util.run_period) == "3") {
        
        if (!_util.use_gpvm)
            flux_file_name = "Systematics/output_rhc_uboone_run0.root";
        else
            flux_file_name = "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/RHC/output_uboone_rhc_run0_merged.root";
        
        
        std::cout << "Using Flux file name: \033[0;31m" << flux_file_name << "\033[0m" <<  std::endl;
        bool boolfile = _util.GetFile(f, flux_file_name );
        if (boolfile == false) gSystem->Exit(0); 
    }
    
    // Get the Flux Histograms
    _util.GetHist(f, h_nue,    "nue/Detsmear/nue_CV_AV_TPC_5MeV_bin");
    _util.GetHist(f, h_nuebar, "nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin");

    double POT_flux = 1.0;

    double POT_Scale_Factor = (flux_scale_factor * _util.config_v.at(_util.k_Run1_Data_POT)) / POT_flux ;

    h_nue->Scale(POT_Scale_Factor);
    h_nuebar->Scale(POT_Scale_Factor);

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    c->SetLeftMargin(0.15);
    c->SetLogy();
    
    // Nue flux
    h_nue->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu}_{e} / cm^{2} / 5 MeV / 0.9 #times 10^{20} POT");
    h_nue->GetXaxis()->CenterTitle();
    h_nue->GetYaxis()->CenterTitle();

    h_nue->GetXaxis()->SetLabelFont(42);
    h_nue->GetXaxis()->SetLabelSize(0.04);
    h_nue->GetYaxis()->SetLabelFont(42);
    h_nue->GetYaxis()->SetLabelSize(0.04);
    h_nue->GetYaxis()->SetTitleSize(0.04);

    h_nue->GetXaxis()->SetRangeUser(0,6);
    // h_nue->GetYaxis()->SetRangeUser(0,150.0e6);
    h_nue->SetLineColor(kBlue+2);
    h_nue->SetFillColor(17);
    
    h_nue->GetXaxis()->SetTitleFont(46);
    h_nue->GetXaxis()->SetTitleSize(18);

    if (!draw_averge) h_nue->Draw("hist");

    TH1D* h_nue_clone = (TH1D*)h_nue->Clone("h_nue_clone");
    
    // Nuebar flux
    h_nuebar->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu}_{e} / cm^{2} / 5 MeV / 0.9 #times 10^{20} POT");
    h_nuebar->GetXaxis()->SetRangeUser(0,6);
    h_nuebar->SetLineColor(kGreen+2);
    h_nuebar->SetFillColor(16);
    if (!draw_averge)h_nuebar->Draw("hist,same");

    // Define the threshold 
    double threshold_energy =  _util.energy_threshold; // Current threshold
    std::cout << "Theshold Energy: " << threshold_energy*1000 << " MeV" << std::endl;
    
    double xbin_th = h_nue->GetXaxis()->FindBin(threshold_energy); // find the x bin to integrate from (threshold)
    double kdar_max = h_nue->GetXaxis()->FindBin( 0.225); // end of KDAR spectrum
    double flux_int_thresh_nue = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1);
    double flux_int_thresh_nuebar = h_nuebar->Integral( xbin_th, h_nuebar->GetNbinsX()+1);

    // Sum the fluxes
    TH1D *summed_flux = (TH1D*)h_nue->Clone("h_summed_flux");

    summed_flux->Add(h_nuebar, 1);

    TH1D* h_summed_flux_clone = (TH1D*)summed_flux->Clone("h_summed_flux_clone");
    
    // Here loop over the cloned histogram of nue, and zero out the flux after the threshold, this is so we can shade out the thresholded area
    for (int p=xbin_th; p < h_summed_flux_clone->GetNbinsX()+1; p++){
        h_summed_flux_clone->SetBinContent(p, 0);
    }

    summed_flux->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e} + #bar{#nu}_{e} / cm^{2} / 5 MeV / 0.9 #times 10^{20} POT");
    summed_flux->GetXaxis()->CenterTitle();
    summed_flux->GetYaxis()->CenterTitle();
    summed_flux->GetXaxis()->SetLabelFont(42);
    summed_flux->GetXaxis()->SetLabelSize(0.04);
    // summed_flux->GetXaxis()->SetTitleSize(0.04);
    summed_flux->GetYaxis()->SetLabelFont(42);
    summed_flux->GetYaxis()->SetLabelSize(0.04);
    summed_flux->GetYaxis()->SetTitleSize(0.04);
    summed_flux->GetXaxis()->SetRangeUser(0,4);
    // summed_flux->GetYaxis()->SetRangeUser(0,170.0e6);
    summed_flux->SetLineColor(kRed+2);
    summed_flux->SetFillColor(17);
    if (draw_averge) summed_flux->Draw("hist");
    

    // Here loop over the cloned histogram of nue, and zero out the flux after the threshold, this is so we can shade out the thresholded area
    for (int p=xbin_th; p < h_nue_clone->GetNbinsX()+1; p++){
        h_nue_clone->SetBinContent(p, 0);
    }

    h_nue_clone->SetLineColor(kBlue+2);
    h_nue_clone->SetLineWidth(0);
    h_nue_clone->SetFillColorAlpha(46, 0.4);
    if (!draw_averge)h_nue_clone->Draw("hist,same");
    // h_nuebar->Draw("hist,same");



    // Get the energy values
    double average_num = 0; // flux numerator
    double average_den = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1) +  h_nuebar->Integral(xbin_th, h_nuebar->GetNbinsX()+1); // flux denominator

    int nue_flux_bins  = h_nue->GetNbinsX();
    int anue_flux_bins = h_nuebar->GetNbinsX();

    for (int k = 1; k < nue_flux_bins+1; k++ ){
        if (k < xbin_th) continue; // skip the thresholded bins
        average_num = average_num + (h_nue->GetBinContent(k) * h_nue->GetBinCenter(k));
    }

    for (int k = 1; k < anue_flux_bins+1; k++ ){
        if (k < xbin_th) continue; // skip the thresholded bins
        average_num = average_num + (h_nuebar->GetBinContent(k) * h_nuebar->GetBinCenter(k));
    }

    double average = average_num / average_den;
    std::cout << "Average Energy: " << average << std::endl;
    

    std::cout << "Total Nue Flux: " << h_nue->Integral() << std::endl;
    std::cout << "Total Nuebar Flux: " << h_nuebar->Integral() << std::endl;
    std::cout << "Total Nue Flux with threshold: " << flux_int_thresh_nue << std::endl;
    std::cout << "Total Nuebar Flux with threshold: " << flux_int_thresh_nuebar << std::endl;

    double sum_flux_thresh_to_KDAR_end = h_nue->Integral( xbin_th, kdar_max) + h_nuebar->Integral( xbin_th, kdar_max);

    h_summed_flux_clone->SetLineWidth(0);
    h_summed_flux_clone->SetFillColorAlpha(46, 0.4);
    h_summed_flux_clone->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e} + #bar{#nu}_{e} / cm^{2} / 0.9 #times 10^{20} POT");
    if (draw_averge) h_summed_flux_clone->Draw("hist,same");

    
    double summed_flux_integral = summed_flux->Integral( xbin_th, summed_flux->GetNbinsX()+1);
    
    std::cout <<  "summed_flux_integral: " << summed_flux_integral << std::endl;
    std::cout << "Fraction of flux from thresh to KDAR end point: " << 100 * sum_flux_thresh_to_KDAR_end / summed_flux_integral << std::endl;

    // Draw the Legend
    TLegend *leg;
    if (!draw_averge)leg = new TLegend(0.27, 0.7, 0.48, 0.9);
    if (draw_averge)leg  = new TLegend(0.27, 0.7, 0.48, 0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (!draw_averge)leg->AddEntry(h_nue,             "NuMI #nu_{e} Flux",                       "f");
    if (!draw_averge)leg->AddEntry(h_nuebar,          "NuMI #bar{#nu}_{e} Flux",                 "f");
    if (draw_averge)leg->AddEntry(summed_flux,        "NuMI #nu_{e} + #bar{#nu}_{e} Flux",       "f");
    if (!draw_averge)leg->Draw();

    // Draw MicroBooNE Simualtion
    _util.Draw_ubooneSim(c, 0.37, 0.92, 0.37, 0.9);

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    if (!draw_averge)c->Print("plots/Integrated_Flux_Separate.pdf");
    if (draw_averge)c->Print("plots/Integrated_Flux_Average.pdf");

}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotTrueVar(){

    // Load in the root file
    TFile * f_mc;
    TTree * mc_tree;      // MC   Tree

    // Get the TTree
    _util.GetFile(f_mc, "../ntuples/neutrinoselection_filt_run1_overlay_weight.root"); // Get the run 1 MC file
    _util.GetTree(f_mc, mc_tree, "nuselection/NeutrinoSelectionFilter");

    // These are the variables we need
    float true_nu_vtx_sce_x, true_nu_vtx_sce_y, true_nu_vtx_sce_z;
    float nu_purity_from_pfp;
    int npi0;
    int ccnc;
    int   nu_pdg;
    float nu_e;
    float elec_e;
    float pi0_e;
    float weightSplineTimesTune;
    float ppfx_cv;    // Weight from PPFX CV
    int   n_showers;
    int   nslice;

    mc_tree->SetBranchAddress("true_nu_vtx_sce_x", &true_nu_vtx_sce_x);
    mc_tree->SetBranchAddress("true_nu_vtx_sce_y", &true_nu_vtx_sce_y);
    mc_tree->SetBranchAddress("true_nu_vtx_sce_z", &true_nu_vtx_sce_z);
    mc_tree->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp);
    mc_tree->SetBranchAddress("npi0", &npi0);
    mc_tree->SetBranchAddress("ccnc",   &ccnc);
    mc_tree->SetBranchAddress("nu_pdg", &nu_pdg);
    mc_tree->SetBranchAddress("nu_e", &nu_e);
    mc_tree->SetBranchAddress("elec_e", &elec_e);
    mc_tree->SetBranchAddress("pi0_e", &pi0_e);
    mc_tree->SetBranchAddress("weightSplineTimesTune",      &weightSplineTimesTune);
    mc_tree->SetBranchAddress("ppfx_cv",                    &ppfx_cv);
    mc_tree->SetBranchAddress("n_showers", &n_showers);
    mc_tree->SetBranchAddress("nslice", &nslice);

    std::vector<std::string> vars = {"nu_e", "elec_e"};

    // Create histograms for the hit purity
    std::vector<std::vector<TH2D*>> h_hit_pur;
    
    // 1D pi0 momentum
    TH1D *h_pi0_momentum = new TH1D("h_true_pi0_momentum", "; #pi^{0} Momentum [GeV/c]; Entries", 40, 0, 2.0);
    
    // 2D shower multiplicity vd nue/electron energy
    TH2D *h_shr_multi_nue_E         = new TH2D("h_shr_multi_nue_E", "; Shower Multiplicty;#nu_{e} Energy [GeV] ", 6, 0, 6, 15, 0, 4.0);
    TH2D *h_shr_multi_elec_e        = new TH2D("h_shr_multi_elec_e", "; Shower Multiplicty;Electron Energy [GeV] ", 6, 0, 6, 15, 0, 4.0);
    TH2D *h_shr_multi_nuebar_E      = new TH2D("h_shr_multi_nuebar_E", "; Shower Multiplicty;#bar{#nu}_{e} Energy [GeV] ", 6, 0, 6, 15, 0, 4.0);
    TH2D *h_shr_multi_elec_e_nuebar = new TH2D("h_shr_multi_elec_e_nuebar", "; Shower Multiplicty;Positron Energy [GeV] ", 6, 0, 6, 15, 0, 4.0);
    
    // Resize hit purity 
    h_hit_pur.resize(vars.size());
    for (unsigned int var = 0; var < h_hit_pur.size(); var++){
        h_hit_pur.at(var).resize(_util.k_classifications_MAX);
    }

    // Create the histograms
    for (unsigned int var = 0; var < h_hit_pur.size(); var++){
        for (unsigned int h = 0; h < h_hit_pur.at(var).size(); h++){
            h_hit_pur.at(var).at(h) = new TH2D(Form("h_hit_pur_%s_%s", vars.at(var).c_str() ,_util.classification_dirs.at(h).c_str()), "", 20, 0, 5.0, 50, 0, 1.1 );
        }
    }

    int mc_tree_total_entries = mc_tree->GetEntries();
    std::cout << "Total MC Events:         " << mc_tree_total_entries << std::endl;

    // Event loop
    for (int ievent = 0; ievent < mc_tree_total_entries; ievent++){

        // See if we want to process all the events
        if (_util.num_events > 0){
            if (ievent >= _util.num_events) break;
        }

        // Alert the user
        if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
    
        // Get the entry in the tree
        mc_tree->GetEntry(ievent); 

        double weight = 1.0;

        // Get the tune weight
        if (_util.weight_tune) weight = weightSplineTimesTune;
        
        // Catch infinate/nan/unreasonably large tune weights
        _util.CheckWeight(weight);

        // Get the PPFX CV flux correction weight
        double weight_flux = 1.0;
        if (_util.weight_ppfx) weight_flux = ppfx_cv;

        _util.CheckWeight(weight_flux);

        if (_util.weight_ppfx) weight = weight * weight_flux;

        // Get the classification
        std::pair<std::string, int> classification = Classify(true_nu_vtx_sce_x, true_nu_vtx_sce_y, true_nu_vtx_sce_z, nu_pdg, ccnc, nu_purity_from_pfp, npi0);      // Classification of the event
        
        // True nue energy
        h_hit_pur.at(0).at(classification.second)->Fill(nu_e, nu_purity_from_pfp, weight);
        
        // True electron energy
        h_hit_pur.at(1).at(classification.second)->Fill(elec_e, nu_purity_from_pfp, weight);

        // Pi0 Momentum
        bool is_in_fv = _util.in_fv(true_nu_vtx_sce_x, true_nu_vtx_sce_y, true_nu_vtx_sce_z);
        if (is_in_fv) h_pi0_momentum->Fill(std::sqrt(pi0_e*pi0_e - 0.134*0.134), weight);

        // Nue cc
        if (nslice == 1 && nu_pdg == 12 && is_in_fv && nu_purity_from_pfp > 0.5 && ccnc == _util.k_CC){

            h_shr_multi_nue_E->Fill(n_showers, nu_e, weight);
            h_shr_multi_elec_e->Fill(n_showers, elec_e, weight);
        }
        // nuebar cc
        if (nslice == 1 && nu_pdg == -12 && is_in_fv && nu_purity_from_pfp > 0.5 && ccnc == _util.k_CC){
            h_shr_multi_nuebar_E->Fill(n_showers, nu_e, weight);
            h_shr_multi_elec_e_nuebar->Fill(n_showers, elec_e, weight);
        }

        

    }

    // Create the resolutions directory for saving the plots to
    _util.CreateDirectory("HitPurity");

    // Now create the directory and save the histograms to file
    // Create the histograms
    for (unsigned int var = 0; var < h_hit_pur.size(); var++){
        for (unsigned int h = 0; h < h_hit_pur.at(var).size(); h++){
            
            TCanvas * c = new TCanvas("c", "c", 500, 500);
            c->SetTopMargin(0.11);

            h_hit_pur.at(var).at(h)->SetStats(kFALSE);

            _util.IncreaseLabelSize(h_hit_pur.at(var).at(h), c);

            gStyle->SetPalette(kBlueGreenYellow);

            if (vars.at(var) == "nu_e") h_hit_pur.at(var).at(h)->SetTitle(Form("%s; True #nu_{e} + #bar{#nu}_{e} Energy [GeV]; Slice Hit Purity", _util.classification_dirs.at(h).c_str()));
            if (vars.at(var) == "elec_e") h_hit_pur.at(var).at(h)->SetTitle(Form("%s; True Electron Energy [GeV]; Slice Hit Purity", _util.classification_dirs.at(h).c_str()));

            h_hit_pur.at(var).at(h)->Draw("colz");

            c->Print(Form("plots/run%s/HitPurity/hit_purity_%s_%s.pdf", _util.run_period, vars.at(var).c_str() ,_util.classification_dirs.at(h).c_str()));
            delete c;

        }
    }

    // Now save the Pi0 Momentum plot
    _util.CreateDirectory("Truth");
    TCanvas * c = new TCanvas("c", "c", 500, 500);
    c->SetTopMargin(0.11);

    h_pi0_momentum->SetStats(kFALSE);

    _util.IncreaseLabelSize(h_pi0_momentum, c);

    h_pi0_momentum->SetLineColor(kAzure - 6);
    h_pi0_momentum->SetLineWidth(2);
    h_pi0_momentum->Draw("hist,E");

    c->Print(Form("plots/run%s/Truth/h_pi0_momentum.pdf", _util.run_period));

    delete c;

    ColumnNorm(h_shr_multi_nue_E);
    ColumnNorm(h_shr_multi_elec_e);
    ColumnNorm(h_shr_multi_nuebar_E);
    ColumnNorm(h_shr_multi_elec_e_nuebar);


    h_shr_multi_nue_E->GetXaxis()->CenterLabels();
    h_shr_multi_elec_e->GetXaxis()->CenterLabels();
    h_shr_multi_nuebar_E->GetXaxis()->CenterLabels();
    h_shr_multi_elec_e_nuebar->GetXaxis()->CenterLabels();


    // Save the 2D shower multiplicity vs nue/elec energy histograms
    Save2DHists(Form("plots/run%s/Truth/h_shower_multiplicity_vs_nue_E.pdf", _util.run_period), h_shr_multi_nue_E);
    Save2DHists(Form("plots/run%s/Truth/h_shower_multiplicity_vs_elec_E.pdf", _util.run_period), h_shr_multi_elec_e);
    Save2DHists(Form("plots/run%s/Truth/h_shower_multiplicity_vs_nuebar_E.pdf", _util.run_period), h_shr_multi_nuebar_E);
    Save2DHists(Form("plots/run%s/Truth/h_shower_multiplicity_vs_elec_E_nuebar.pdf", _util.run_period), h_shr_multi_elec_e_nuebar);
    

}
// -----------------------------------------------------------------------------
std::pair<std::string, int> UtilityPlotter::Classify(float true_nu_vtx_sce_x, float true_nu_vtx_sce_y, float true_nu_vtx_sce_z, int nu_pdg, int ccnc, float nu_purity_from_pfp, int npi0){
   
    bool is_in_fv = _util.in_fv(true_nu_vtx_sce_x, true_nu_vtx_sce_y, true_nu_vtx_sce_z);

    // Out of Fiducial Volume Event
    if (!is_in_fv) {
        // std::cout << "Purity of out of FV event: "<< nu_purity_from_pfp << std::endl;
        if (nu_purity_from_pfp < 0.0) return std::make_pair("unmatched",_util.k_unmatched);
        else return std::make_pair("nu_out_fv",_util.k_nu_out_fv);
    }
    // In FV event
    else {

        // Charged Current 
        if (ccnc == _util.k_CC){

            // NuMu CC
            if (nu_pdg == 14 || nu_pdg == -14){

                // Purity is low so return cosmic
                if (nu_purity_from_pfp < 0.0)return std::make_pair("unmatched",_util.k_unmatched);
                
                if (npi0 > 0) return std::make_pair("numu_cc_pi0", _util.k_numu_cc_pi0); // has a pi0
                else return std::make_pair("numu_cc",_util.k_numu_cc);

            }
            // Nue CC
            else if (nu_pdg == 12){
                
                if (nu_purity_from_pfp > 0.0)                                 return std::make_pair("nue_cc",       _util.k_nue_cc);    // purity > 0.5% so signal
                else                                                          return std::make_pair("unmatched_nue",_util.k_unmatched_nue); // These events were not picked up by pandora at all

            }
            else if (nu_pdg == -12){
                
                if (nu_purity_from_pfp > 0.0)                                  return std::make_pair("nuebar_cc",       _util.k_nuebar_cc); // purity > 0.5% so signal
                else                                                           return std::make_pair("unmatched_nuebar",_util.k_unmatched_nuebar); // These events were not picked up by pandora at all

            }
            // Unknown Neutrino Type
            else {
                std::cout << "Unknown Neutrino Type..., This will also mess up the efficecy if this occurs!" << std::endl;
                return std::make_pair("unmatched",_util.k_unmatched);
            }

        }
        // Neutral Current
        else {

            // Purity is low so return cosmic
            if (nu_purity_from_pfp < 0) return std::make_pair("unmatched",_util.k_unmatched);

            if (npi0 > 0) return std::make_pair("nc_pi0",_util.k_nc_pi0);
            else return std::make_pair("nc",_util.k_nc);
        }
    
    } // End if in FV


}
// -----------------------------------------------------------------------------
void UtilityPlotter::Save2DHists(const char* printname, TH2D* hist){

    TCanvas * c = new TCanvas("c", "c", 500, 500);
    c->SetTopMargin(0.11);

    hist->SetStats(kFALSE);

    _util.IncreaseLabelSize(hist, c);

    gStyle->SetPalette(kBlueGreenYellow);

    hist->Draw("colz");

    c->Print(printname);
    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::ColumnNorm(TH2D* hist){

    // Loop over rows
    for (int row=1; row<hist->GetXaxis()->GetNbins()+1; row++) {
        double integral = 0;

        // Loop over columns and get the integral
        for (int col=1; col<hist->GetYaxis()->GetNbins()+1; col++){
            integral+=hist->GetBinContent(row, col);            
        }

        // Now normalise the column entries by the integral
        for (int col=1; col<hist->GetYaxis()->GetNbins()+1; col++){
            hist->SetBinContent(row,col, hist->GetBinContent(row, col)/ integral );
            
        }
    } 

}
// -----------------------------------------------------------------------------
void UtilityPlotter::StudyPPFXWeights(){

    std::vector<TH1D*> v_hist;
    v_hist.resize(_util.k_classifications_MAX);

    for (unsigned int c = 0; c < v_hist.size(); c++){
        v_hist.at(c) = new TH1D(Form("h_ppfx_weight_%s", _util.classification_dirs.at(c).c_str()), "; PPFX Weights All Universes; Entries", 100, 0, 2.0);
    }


    // Loop over the entries in the TTree
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 

        std::vector<double> vec_universes;

        // Convert from unsigned short to double and push back -- divide by 1000 to undo previous *1000
        for (unsigned int j = 0; j < weightsPPFX->size(); j++){
            vec_universes.push_back( (double) weightsPPFX->at(j)/1000.0);
        }

        // Now loop over the universes
        for (unsigned int uni = 0; uni < vec_universes.size(); uni++){

            // Background events
            if ( *classification == "nu_out_fv"){
                v_hist.at(_util.k_nu_out_fv)->Fill(vec_universes.at(uni));
            }
            if ( *classification == "cosmic"){
                v_hist.at(_util.k_cosmic)->Fill(vec_universes.at(uni));
            }
            if ( *classification == "numu_cc"){
                v_hist.at(_util.k_numu_cc)->Fill(vec_universes.at(uni));
            }
            if ( *classification == "numu_cc_pi0"){
                v_hist.at(_util.k_numu_cc_pi0)->Fill(vec_universes.at(uni));
            }
            if ( *classification == "nc"){
                v_hist.at(_util.k_nc)->Fill(vec_universes.at(uni));
            }
            if ( *classification == "nc_pi0"){
                v_hist.at(_util.k_nc_pi0)->Fill(vec_universes.at(uni));
            }
        }
            
        // Generated event
        // if ( (*classification == "nue_cc"|| *classification == "nuebar_cc" || *classification == "unmatched_nue" || *classification == "cosmic_nue" || *classification == "unmatched_nuebar" || *classification == "cosmic_nuebar") && gen == true) {}
        // Signal event
        // if ((*classification == "nue_cc" || *classification == "nuebar_cc" || *classification == "unmatched_nue" || *classification == "unmatched_nuebar") && gen == false) {}
     
    }

    TCanvas * c = new TCanvas("c", "c", 500, 500);
    c->SetTopMargin(0.11);

     TLegend *leg = new TLegend(0.2, 0.55, 0.4, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (unsigned int h = 0; h < v_hist.size(); h++){

        if (h == _util.k_nue_cc || 
            h == _util.k_nuebar_cc || 
            h == _util.k_unmatched || 
            h == _util.k_unmatched_nue || 
            h == _util.k_cosmic_nue || 
            h == _util.k_unmatched_nuebar || 
            h == _util.k_cosmic_nuebar || 
            h == _util.k_leg_ext || 
            h == _util.k_leg_data || 
            h == _util.k_leg_dirt) {
                continue;
        }

        v_hist.at(h)->Scale(1.0/v_hist.at(h)->Integral());

        v_hist.at(h)->SetStats(kFALSE);

        _util.IncreaseLabelSize(v_hist.at(h), c);
        v_hist.at(h)->SetMaximum(0.06);

        v_hist.at(h)->GetYaxis()->SetTitleOffset(1.5);

        leg->AddEntry(v_hist.at(h),Form("%s", _util.classification_dirs.at(h).c_str()), "l");

        if ( h == _util.k_nu_out_fv ){
            v_hist.at(h)->SetLineColor(30);
        }
        if ( h == _util.k_cosmic     ){
            v_hist.at(h)->SetLineColor(38);
        }
        if ( h == _util.k_numu_cc    ){
            v_hist.at(h)->SetLineColor(28);
        }
        if ( h == _util.k_numu_cc_pi0){
            v_hist.at(h)->SetLineColor(4);
        }
        if ( h == _util.k_nc         ){
            v_hist.at(h)->SetLineColor(36);
        }
        if ( h ==  _util.k_nc_pi0){
            v_hist.at(h)->SetLineColor(1);
        }

        v_hist.at(h)->SetLineWidth(2);
        v_hist.at(h)->Draw("hist,same");

    }
    leg->Draw();
    
    // Create the Bins directory for saving the plots to
    _util.CreateDirectory("Systematics/Misc");
    c->Print("plots/run1/Systematics/Misc/ppfx_weights_backgrounds.pdf");

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareEfficiency(){

    TH1D* h_eff_run1 = new TH1D("h_efficiency_run1", "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);
    TH1D* h_pur_run1 = new TH1D("h_purity_run1",     "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);
    TH1D* h_eff_clone_run1 = new TH1D("h_eff_err_run1",     "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);;

    TH1D* h_eff_run3 = new TH1D("h_efficiency_run3", "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);
    TH1D* h_pur_run3 = new TH1D("h_purity_run3",     "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);
    TH1D* h_eff_clone_run3 = new TH1D("h_eff_err_run3",     "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);;

    PopulateEff(h_eff_run1, h_pur_run1, h_eff_clone_run1,  "files/trees/nuexsec_selected_tree_mc_run1.root");
    PopulateEff(h_eff_run3, h_pur_run3, h_eff_clone_run3, "files/trees/nuexsec_selected_tree_mc_run3.root");
    
    
    TCanvas *c = new TCanvas("c", "c", 600, 500);
    c->SetGridy();
    c->SetBottomMargin(0.12);

    TLegend *leg_stack = new TLegend(0.59, 0.89, 0.89, 0.69);
    leg_stack->SetBorderSize(0);
    // leg_stack->SetFillStyle(0);
    leg_stack->AddEntry(h_eff_run1, "Efficiency Run1","ELP");
    leg_stack->AddEntry(h_eff_run3, "Efficiency Run3","ELP");
    leg_stack->AddEntry(h_pur_run1, "Purity Run1",    "lp");
    leg_stack->AddEntry(h_pur_run3, "Purity Run3",    "lp");

    h_eff_run1->Draw("LP");
    h_pur_run1->Draw("LP,same");
    

    h_eff_run3->SetLineColor(kAzure+5);
    h_eff_run3->Draw("LP,same");

    h_eff_clone_run3->SetLineColor(kAzure+5);

    h_pur_run3->SetLineColor(kViolet-5);
    h_pur_run3->Draw("LP,same");
    
    // Draw vertical lines to help the eye
    TLine *line;
    for (unsigned int l=1; l < _util.k_cuts_MAX+1; l++){
        line  = new TLine( h_eff_run1->GetBinCenter(l) ,   0 , h_eff_run1->GetBinCenter(l)  ,  1.1);
        line->SetLineColor(12);
        line->SetLineStyle(kDotted);
        line->Draw();
    }

    leg_stack->Draw();

    h_eff_clone_run1->Draw("E,X0,same");
    h_eff_clone_run3->Draw("E,X0, same");

    h_eff_run1->GetXaxis()->SetTickLength(0.00);
    h_pur_run1->GetXaxis()->SetTickLength(0.00);

    c->Print("plots/efficiency_run1_run3_comparison.pdf");


    // Now lets plot the difference between the two runs

    TH1D* h_eff_diff = (TH1D*) h_eff_run1->Clone("h_eff_diff");
    h_eff_diff->Add(h_eff_run3, -1);
    
    // Set the bin error by adding in quadrature
    for (int bin = 1; bin < h_eff_diff->GetNbinsX()+1; bin++){
        double run1_err = h_eff_clone_run1->GetBinError(bin);
        double run3_err = h_eff_clone_run3->GetBinError(bin);
        double tot_err = std::sqrt(run1_err*run1_err + run3_err*run3_err);
        h_eff_diff->SetBinError(bin, tot_err);
    }

    h_eff_diff->SetStats(kFALSE);
    h_eff_diff->SetMarkerStyle(20);
    h_eff_diff->SetMarkerSize(0.5);
    h_eff_diff->SetLineWidth(2);
    h_eff_diff->GetYaxis()->SetTitle("Difference \%");
    h_eff_diff->Scale(100);
    h_eff_diff->GetYaxis()->SetRangeUser(-8, 8);
    h_eff_diff->GetXaxis()->SetLabelFont(62);
    h_eff_diff->GetXaxis()->SetLabelSize(0.03);

    TH1D* h_pur_diff = (TH1D*) h_pur_run1->Clone("h_pur_diff");
    h_pur_diff->Add(h_pur_run3, -1);
    h_pur_diff->SetLineColor(kRed+2);
    h_pur_diff->SetStats(kFALSE);
    h_pur_diff->SetMarkerStyle(20);
    h_pur_diff->SetMarkerSize(0.5);
    h_pur_diff->SetLineWidth(2);
    h_pur_diff->Scale(100);

    TH1D* h_eff_diff_clone = (TH1D*) h_eff_diff->Clone();

    for (int bin = 1; bin < h_eff_diff->GetNbinsX()+1; bin++){
        h_eff_diff->SetBinError(bin, 0.0);
    }


    TCanvas *c2 = new TCanvas("c2", "c2", 600, 500);
    c2->SetGridy();
    c2->SetBottomMargin(0.12);
    h_eff_diff->Draw("LP");
    h_eff_diff_clone->Draw("E,X0,same");
    h_pur_diff->Draw("LP, same");

    for (unsigned int l=1; l < _util.k_cuts_MAX+1; l++){
        line  = new TLine( h_eff_run1->GetBinCenter(l) ,   -8 , h_eff_run1->GetBinCenter(l)  ,  8);
        line->SetLineColor(12);
        line->SetLineStyle(kDotted);
        line->Draw();
    }

    TLegend *leg_stack2 = new TLegend(0.2, 0.89, 0.4, 0.69);
    leg_stack2->SetBorderSize(0);
    // leg_stack2->SetFillStyle(0);
    leg_stack2->AddEntry(h_eff_diff, "Efficiency","ELP");
    leg_stack2->AddEntry(h_pur_diff, "Purity","lp");
    leg_stack2->Draw();

    c2->Print("plots/efficiency_run1_run3_comparison_diff.pdf");


    // Now lets make a relative differnce plot
    TH1D* h_eff_diff_rel = new TH1D("h_efficiency_diff_rel", "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);
    TH1D* h_pur_diff_rel = new TH1D("h_purity_diff_rel",     "", _util.k_cuts_MAX, 0, _util.k_cuts_MAX);

    for (int k=0; k < h_eff_run1->GetNbinsX();k++){
        if ( k == 0 ){
            h_eff_diff_rel ->Fill(_util.cut_dirs_pretty.at(k).c_str(), 0);
            h_pur_diff_rel ->Fill(_util.cut_dirs_pretty.at(k).c_str(), 0);
        }
        else {
            h_eff_diff_rel ->Fill(_util.cut_dirs_pretty.at(k).c_str(), h_eff_diff->GetBinContent(k+1) - h_eff_diff->GetBinContent(k));
            h_pur_diff_rel ->Fill(_util.cut_dirs_pretty.at(k).c_str(), h_pur_diff->GetBinContent(k+1) - h_pur_diff->GetBinContent(k));
        }
        
        // Now set the bin error
        double err1 = h_eff_diff_clone->GetBinError(k+1);
        double err2 = h_eff_diff_clone->GetBinError(k);
        double tot_err = std::sqrt(err1*err1 + err2*err2);
        h_eff_diff_rel->SetBinError(k, tot_err);

        h_pur_diff_rel->SetBinError(k+1, 0);
    }

    h_eff_diff_rel->SetStats(kFALSE);
    h_eff_diff_rel->SetMarkerStyle(20);
    h_eff_diff_rel->SetMarkerSize(0.5);
    h_eff_diff_rel->SetLineWidth(2);
    h_eff_diff_rel->GetYaxis()->SetTitle("Rel. Difference (Run 1 - Run 3) \%");
    // h_eff_diff_rel->Scale(100);
    h_eff_diff_rel->GetYaxis()->SetRangeUser(-3, 3);

    h_eff_diff_rel->GetXaxis()->SetLabelFont(62);
    h_eff_diff_rel->GetXaxis()->SetLabelSize(0.03);

    h_pur_diff_rel->SetLineColor(kRed+2);
    h_pur_diff_rel->SetStats(kFALSE);
    h_pur_diff_rel->SetMarkerStyle(20);
    h_pur_diff_rel->SetMarkerSize(0.5);
    h_pur_diff_rel->SetLineWidth(2);
    // h_pur_diff_rel->Scale(100);

    TH1D* h_eff_diff_rel_clone = (TH1D*) h_eff_diff_rel->Clone();

    for (int bin = 1; bin < h_eff_diff->GetNbinsX()+1; bin++){
        h_eff_diff_rel->SetBinError(bin, 0.0);
    }

    TCanvas *c3 = new TCanvas("c3", "c3", 600, 500);
    c3->SetGridy();
    c3->SetBottomMargin(0.12);
    h_eff_diff_rel->Draw("LP");
    h_pur_diff_rel->Draw("LP, same");
    h_eff_diff_rel_clone->Draw("E, X0, same");

    // Draw vertical lines to help the eye
    for (unsigned int l=1; l < _util.k_cuts_MAX+1; l++){
        line  = new TLine( h_eff_run1->GetBinCenter(l) ,   -3 , h_eff_run1->GetBinCenter(l)  ,  3);
        line->SetLineColor(12);
        line->SetLineStyle(kDotted);
        line->Draw();
    }

    leg_stack2->Draw();

    c3->Print("plots/efficiency_run1_run3_comparison_reldiff.pdf");



}
// -----------------------------------------------------------------------------
void UtilityPlotter::PopulateEff(TH1D* h_eff, TH1D *h_pur, TH1D* h_eff_clone, const char* input_file){

    TTree *mc_tree;

    TFile *f;

    // The uglyest hardcoded monstosity known to the universe. SORT this out KRISH...
    f = TFile::Open(input_file);

    std::vector<double> efficiency_v; // efficiency vector
    std::vector<double> eff_err_v; // efficiency error vector
    std::vector<double> purity_v    ; // purity vector

    double efficiency, purity, eff_err;

    _util.GetTree(f, mc_tree, "mc_eff_tree");
    mc_tree->SetBranchAddress("efficiency", &efficiency);
    mc_tree->SetBranchAddress("purity",     &purity);
    mc_tree->SetBranchAddress("eff_err",    &eff_err);

    int num_entries = mc_tree->GetEntries();

    // Fill the efficiency and purity vectors
    for (int y=0; y < num_entries; y++){
        mc_tree->GetEntry(y); 
        efficiency_v.push_back(efficiency);
        purity_v.push_back(purity);
        eff_err_v.push_back(eff_err);
    }

    for (unsigned int k=0; k < efficiency_v.size();k++){
        if (k == 0 || k == 1) {
            h_eff ->Fill(_util.cut_dirs_pretty.at(k).c_str(), 1.0); // We need to put these back to 1 (becasue of the way we define things after slice id)
            h_eff_clone->Fill(_util.cut_dirs_pretty.at(k).c_str(), 1.0);
        }
        else {
            h_eff ->Fill(_util.cut_dirs_pretty.at(k).c_str(), efficiency_v.at(k));
            h_eff_clone->Fill(_util.cut_dirs_pretty.at(k).c_str(), efficiency_v.at(k));
        }
        
        h_pur ->Fill(_util.cut_dirs_pretty.at(k).c_str(), purity_v.at(k));
        
        h_eff->SetBinError(k+1, 0);
        h_pur->SetBinError(k+1, 0);
        h_eff_clone->SetBinError(k+1, eff_err_v.at(k));
    }

    h_eff->GetYaxis()->SetRangeUser(0, 1.1);
    h_eff->SetStats(kFALSE);
    h_eff->SetMarkerStyle(20);
    h_eff->SetMarkerSize(0.5);
    h_eff->SetLineWidth(2);
    h_eff->GetXaxis()->SetTitleSize(0.05);
    h_eff->GetYaxis()->SetLabelSize(0.05);
    h_eff->GetYaxis()->SetTitleSize(0.05);
    h_eff->GetXaxis()->SetLabelFont(62);
    h_eff->GetXaxis()->SetLabelSize(0.03);
    h_eff_clone->SetLineWidth(2);

    h_pur->SetLineColor(kRed+2);
    h_pur->SetStats(kFALSE);
    h_pur->SetMarkerStyle(20);
    h_pur->SetMarkerSize(0.5);
    h_pur->SetLineWidth(2);


}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareDetVarEfficiency(){

    TFile *f_cv, *f_cv_intrinsic;

    gStyle->SetOptStat(0);


    // Get the CV files
    f_cv = TFile::Open( Form("files/nuexsec_mc_run%s_CV.root", _util.run_period ));
    f_cv_intrinsic = TFile::Open( Form("files/nuexsec_mc_run%s_CV_intrinsic.root", _util.run_period ));

    // Now get the efficiency denominators and numerators
    TH1D *h_cv_den, *h_cv_intrinsic_den;
    TH1D *h_cv_num, *h_cv_intrinsic_num;

    _util.GetHist(f_cv, h_cv_den, "TEff/h_true_elec_E_rebin_Unselected");
    _util.GetHist(f_cv, h_cv_num, "TEff/h_true_elec_E_rebin_dEdx_max_no_tracks");

    _util.GetHist(f_cv_intrinsic, h_cv_intrinsic_den, "TEff/h_true_elec_E_rebin_Unselected");
    _util.GetHist(f_cv_intrinsic, h_cv_intrinsic_num, "TEff/h_true_elec_E_rebin_dEdx_max_no_tracks");


    std::vector<double> err_cv(h_cv_den->GetNbinsX(), 0.0);
    std::vector<double> err_cv_intrinsic(h_cv_den->GetNbinsX(), 0.0);

    // Get the errors for using binomial errors
    for (int bin = 1; bin < h_cv_den->GetNbinsX()+1; bin ++){
        double n = h_cv_num->GetBinContent(bin);
        double N = h_cv_den->GetBinContent(bin);
        err_cv.at(bin-1) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

        

        n = h_cv_intrinsic_num->GetBinContent(bin)*(1.0/0.124933);
        N = h_cv_intrinsic_den->GetBinContent(bin)*(1.0/0.124933);
        err_cv_intrinsic.at(bin-1) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

        std::cout << N << " " << n<< " "<< err_cv_intrinsic.at(bin-1) << std::endl;
    
    }
    


    // Now divide the histograms
    h_cv_num->Divide(h_cv_den);
    h_cv_intrinsic_num->Divide(h_cv_intrinsic_den);


    // Now set the bin error
    for (int bin = 1; bin < h_cv_den->GetNbinsX()+1; bin ++){
        h_cv_num->SetBinError(bin, err_cv.at(bin-1));
        h_cv_intrinsic_num->SetBinError(bin, err_cv_intrinsic.at(bin-1));
    }

    // And save them
    TCanvas *c = new TCanvas("c", "c", 500, 500);

    TLegend *leg = new TLegend(0.15, 0.7, 0.55, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    h_cv_num->SetTitle(";True Electron Energy [GeV]; Efficiency");
    h_cv_num->GetYaxis()->SetTitleOffset(1.4);
    h_cv_num->SetMaximum(0.5);
    h_cv_num->SetMinimum(0);
    h_cv_num->SetLineWidth(2);
    h_cv_num->SetLineColor(kBlack);
    h_cv_intrinsic_num->SetLineWidth(2);
    
    h_cv_num->Draw("hist,E");
    h_cv_intrinsic_num->Draw("hist,E,same");

    leg->AddEntry(h_cv_num, "DetVar CV", "l");
    leg->AddEntry(h_cv_intrinsic_num, "Intrinsic DetVar CV", "l");
    leg->Draw();


    c->Print("plots/run1/detvar/comparisons/Eff_CV_comparison.pdf");

    delete c;


}
// -----------------------------------------------------------------------------
void UtilityPlotter::TestModelDependence(){

    gStyle->SetOptStat(0);

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH2D* h_temp_2D;
    TH1D* h_temp;

    // The covariance matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/er/h_cov_tot_mcxsec_reco",_util.xsec_var));
    TH2D* h_cov = (TH2D*)h_temp_2D->Clone();
    h_cov->SetDirectory(0);

    // Data XSec
    h_temp  = (TH1D*)fxsec->Get(Form("%s/er/h_data_xsec_stat_sys_reco", _util.xsec_var));
    TH1D* h_dataxsec = (TH1D*)h_temp->Clone();
    h_dataxsec->SetDirectory(0);
    h_dataxsec->SetLineColor(kBlack);
    h_dataxsec->SetMarkerStyle(20);
    h_dataxsec->SetMarkerSize(0.5);

    h_temp  = (TH1D*)fxsec->Get(Form("%s/er/h_data_xsec_stat_reco", _util.xsec_var));
    TH1D* h_dataxsec_stat = (TH1D*)h_temp->Clone();
    h_dataxsec_stat->SetDirectory(0);
    h_dataxsec_stat->SetLineColor(kBlack);
    h_dataxsec_stat->SetMarkerStyle(20);
    h_dataxsec_stat->SetMarkerSize(0.5);

    fxsec->Close();

    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Create a vector for the models
    std::vector<std::string> models = {
        "CV",
        "mec",
        "nogtune",
        "nopi0tune",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_CV,
        k_model_mec,
        k_model_nogtune,
        k_model_nopi0tune,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };

    // Create the vector of histograms
    std::vector<TH2D*> h_response_model(models.size());
    std::vector<TH1D*> h_mcxsec_true_model(models.size());
    std::vector<TH1D*> h_mcxsec_reco_model(models.size());


    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
        
        // Response Matrix
        h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run1_%s_0_smearing", models.at(m).c_str(), vars.at(k_var_trueX).c_str(), models.at(m).c_str()));
        if (h_temp_2D == NULL) std::cout <<"Help!" << m << std::endl;
        h_response_model.at(m) = (TH2D*)h_temp_2D->Clone();

        // MC xsec in True
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", models.at(m).c_str(), vars.at(k_var_trueX).c_str(), vars.at(k_var_trueX).c_str()));
        h_mcxsec_true_model.at(m) = (TH1D*)h_temp->Clone();
        h_mcxsec_reco_model.at(m) = (TH1D*)h_dataxsec->Clone();
       
        _util.MatrixMultiply(h_mcxsec_true_model.at(m), h_mcxsec_reco_model.at(m), h_response_model.at(k_model_CV), "true_reco", true);
    }

    // Set the line colours
    h_mcxsec_reco_model.at(k_model_CV)       ->SetLineColor(kRed+2);
    h_mcxsec_reco_model.at(k_model_mec)      ->SetLineColor(kGreen+2);
    h_mcxsec_reco_model.at(k_model_nogtune)  ->SetLineColor(kBlue+2);
    h_mcxsec_reco_model.at(k_model_nopi0tune)->SetLineColor(kPink+1);
    h_mcxsec_reco_model.at(k_model_FLUGG)    ->SetLineColor(kYellow+2);
    h_mcxsec_reco_model.at(k_model_tune1)    ->SetLineColor(kOrange-1);
    

    // Now lets plot
    TCanvas *c = new TCanvas("c", "c", 500, 500);
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.15);
    h_dataxsec->GetYaxis()->SetTitleOffset(1.7);
    // h_mcxsec_reco->SetMaximum(1.5);
    h_dataxsec->Draw("E1,X0,same");

    if (_util.zoom && std::string(_util.xsec_var) == "elec_cang")
        h_dataxsec->GetXaxis()->SetRangeUser(0.6, 1.0);
    
    h_mcxsec_reco_model.at(k_model_CV)->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_mec)->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_nogtune)->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_nopi0tune)->Draw("hist,same");
    // h_mcxsec_reco_model.at(k_model_FLUGG)->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_tune1)->Draw("hist,same");
    h_dataxsec->Draw("E1,X0,same");
    h_dataxsec_stat->Draw("E1,X0,same");

    TLegend *leg = new TLegend(0.4, 0.6, 0.75, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_dataxsec, "Data (Stat. + Sys.)", "ep");
    
    double chi, pval;
    int ndof;
    std::cout << "CV" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_CV), h_dataxsec, h_cov, chi, ndof, pval);
    // _util.CalcChiSquaredNoCorr(h_mcxsec_reco_model.at(k_model_CV), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_CV),  Form("MC (CV) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");
    
    std::cout << "1.5 MEC" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_mec), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_mec),  Form("MC (1.5 #times MEC) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");
    
    std::cout << "No gTune" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_nogtune), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_nogtune),  Form("MC (no gTune) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");
    
    std::cout << "no pi0 tune" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_nopi0tune), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_nopi0tune),  Form("MC (no #pi^{0} Tune) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");
    
    
    // _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_FLUGG), h_dataxsec, h_cov, chi, ndof, pval);
    // leg->AddEntry(h_mcxsec_reco_model.at(k_model_FLUGG),  Form("MC (FLUGG Flux) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");

    std::cout << "Tune 1" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_tune1), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_tune1),  Form("MC (Tune 1) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");
    
    leg->Draw();

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.52, 0.92, 0.52, 0.92);

    if (_util.zoom)
        c->Print(Form("plots/run%s/Models/%s/DataModelComparison_zoom.pdf", _util.run_period, _util.xsec_var));
    else
        c->Print(Form("plots/run%s/Models/%s/DataModelComparison.pdf", _util.run_period, _util.xsec_var));

    fxsec->Close();

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareDataCrossSections(){

    gStyle->SetOptStat(0);

    /// Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH1D* h_temp;

    // Data XSec
    h_temp  = (TH1D*)fxsec->Get(Form("%s/er/h_data_xsec_stat_sys_reco", _util.xsec_var));
    TH1D* h_dataxsec = (TH1D*)h_temp->Clone();
    h_dataxsec->SetDirectory(0);
    h_dataxsec->SetLineColor(kBlack);
    h_dataxsec->SetMarkerStyle(20);
    h_dataxsec->SetMarkerSize(0.5);

    h_temp  = (TH1D*)fxsec->Get(Form("%s/er/h_data_xsec_stat_reco", _util.xsec_var));
    TH1D* h_dataxsec_stat = (TH1D*)h_temp->Clone();
    h_dataxsec_stat->SetDirectory(0);
    h_dataxsec_stat->SetLineColor(kBlack);
    h_dataxsec_stat->SetMarkerStyle(20);
    h_dataxsec_stat->SetMarkerSize(0.5);

    fxsec->Close();

    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");
    
    // Data X Sec MEC
    h_temp  = (TH1D*)fxsec->Get(Form("mec/%s/h_run1_CV_0_%s_data_xsec", vars.at(k_var_recoX).c_str(), vars.at(k_var_recoX).c_str()));
    TH1D* h_datasec_reco_mec = (TH1D*)h_temp->Clone();
    h_datasec_reco_mec->SetLineColor(kGreen+2);
    h_datasec_reco_mec->Scale(1.0, "width");
    h_datasec_reco_mec->SetLineWidth(2);

    // Data X Sec No Genie Tune
    h_temp  = (TH1D*)fxsec->Get(Form("nogtune/%s/h_run1_CV_0_%s_data_xsec", vars.at(k_var_recoX).c_str(), vars.at(k_var_recoX).c_str()));
    TH1D* h_datasec_reco_nogtune = (TH1D*)h_temp->Clone();
    h_datasec_reco_nogtune->SetLineColor(kBlue+2);
    h_datasec_reco_nogtune->Scale(1.0, "width");
    h_datasec_reco_nogtune->SetLineWidth(2);
    
    // Data X Sec Tune 1
    h_temp  = (TH1D*)fxsec->Get(Form("tune1/%s/h_run1_CV_0_%s_data_xsec", vars.at(k_var_recoX).c_str(), vars.at(k_var_recoX).c_str()));
    TH1D* h_datasec_reco_tune1 = (TH1D*)h_temp->Clone();
    h_datasec_reco_tune1->SetLineColor(kOrange-1);
    h_datasec_reco_tune1->Scale(1.0, "width");
    h_datasec_reco_tune1->SetLineWidth(2);

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.15);
    h_dataxsec->GetYaxis()->SetTitleOffset(1.7);
    h_dataxsec->Draw("E1,X0");
    h_datasec_reco_mec->Draw("hist,same");
    h_datasec_reco_nogtune->Draw("hist,same");
    h_datasec_reco_tune1->Draw("hist,same");
    h_dataxsec->Draw("E1,same,X0");
    h_dataxsec_stat->Draw("E1,same,X0");
    

    TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_dataxsec, "Data (Stat. + Sys.)", "ep");
    leg->AddEntry(h_datasec_reco_mec, "Data with 1.5 #times MEC Model", "l");
    leg->AddEntry(h_datasec_reco_nogtune, "Data with no gTune Model", "l");
    leg->AddEntry(h_datasec_reco_tune1, "Data with Tune 1 Model", "l");
    leg->Draw();

    c->Print(Form("plots/run%s/Models/%s/ModelDataComparison.pdf", _util.run_period, _util.xsec_var));

    fxsec->Close();

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareSmearing(){

    gStyle->SetOptStat(0);

    /// Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH1D* h_temp;
    TH2D* h_temp_2D;

    // Get the MC covariance Matrix
    h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/er/h_cov_tot_mcxsec_smear_true", _util.xsec_var ));
    TH2D* h_cov_smear_tot = (TH2D*)h_temp_2D->Clone();
    h_cov_smear_tot->SetDirectory(0);

    // Get the reco xsec
    h_temp = (TH1D*)fxsec->Get(Form("%s/er/h_mc_xsec_reco",_util.xsec_var));
    TH1D* h_mcxsec_reco = (TH1D*)h_temp->Clone();
    h_mcxsec_reco->SetDirectory(0);
    h_mcxsec_reco->SetLineColor(kBlack);


    fxsec->Close();

    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Create a vector for the models
    std::vector<std::string> models = {
        "mec",
        "nogtune",
        "nopi0tune",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_mec,
        k_model_nogtune,
        k_model_nopi0tune,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };

    // Create the vector of histograms
    std::vector<TH2D*> h_response_model(models.size());
    std::vector<TH1D*> h_mcxsec_reco_model(models.size());
    
    // MC Xsec True
    h_temp  = (TH1D*)fxsec->Get(Form("CV/%s/h_run1_CV_0_%s_mc_xsec",vars.at(k_var_trueX).c_str(), vars.at(k_var_trueX).c_str()));
    TH1D* h_mcxsec_true         = (TH1D*)h_temp->Clone();
    
    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
        
        // clone to get the binning, these hists get reset
        h_mcxsec_reco_model.at(m) = (TH1D*)h_mcxsec_reco->Clone();

        // Get the response matrix
        h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run1_%s_0_smearing", models.at(m).c_str(), vars.at(k_var_trueX).c_str(), models.at(m).c_str()));
        h_response_model.at(m) = (TH2D*)h_temp_2D->Clone();

        // Apply the response matrix to the CV MC True dist
        _util.MatrixMultiply(h_mcxsec_true, h_mcxsec_reco_model.at(m), h_response_model.at(m), "true_reco", true);
    }

    // Set the line colours
    h_mcxsec_reco_model.at(k_model_mec)      ->SetLineColor(kGreen+2);
    h_mcxsec_reco_model.at(k_model_nogtune)  ->SetLineColor(kBlue+2);
    h_mcxsec_reco_model.at(k_model_nopi0tune)->SetLineColor(kPink+1);
    h_mcxsec_reco_model.at(k_model_FLUGG)    ->SetLineColor(kYellow+2);
    h_mcxsec_reco_model.at(k_model_tune1)    ->SetLineColor(kOrange-1);

    // Set the Bin errors for the MC truth
    for (int bin = 1; bin < h_mcxsec_reco->GetNbinsX()+1; bin++){
        
        // Set the bin error of the CV to be the stat plus smearing uncertainty 
        double err_mc_true = (h_mcxsec_true->GetBinError(bin) / h_mcxsec_true->GetBinContent(bin)) * h_mcxsec_reco->GetBinContent(bin);
        h_mcxsec_reco->SetBinError(bin, err_mc_true + std::sqrt(h_cov_smear_tot->GetBinContent(bin, bin)));        
        
        // Set the bin error for each model to be ~2%
        for (unsigned int m = 0; m < models.size(); m++){
            h_mcxsec_reco_model.at(m)->SetBinError(bin, 0.02*h_mcxsec_reco_model.at(m)->GetBinContent(bin));
        }
    }

    if (_util.zoom && std::string(_util.xsec_var) == "elec_cang")
        h_mcxsec_reco->GetXaxis()->SetRangeUser(0.6, 1.0);

    // Now lets plot
    TCanvas *c = new TCanvas("c", "c", 500, 500);
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.15);
    h_mcxsec_reco->GetYaxis()->SetTitleOffset(1.7);
    h_mcxsec_reco->Draw("hist,E");
    
    // Choose what models to draw
    h_mcxsec_reco_model.at(k_model_mec)      ->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_nogtune)  ->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_nopi0tune)->Draw("hist,same");
    // h_mcxsec_reco_model.at(k_model_FLUGG)    ->Draw("hist,E,same");
    h_mcxsec_reco_model.at(k_model_tune1)    ->Draw("hist,E,same");
    h_mcxsec_reco->Draw("hist,E,same");

    // Create the legend
    TLegend *leg = new TLegend(0.4, 0.5, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_mcxsec_reco, "MC (Stat.)", "el");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_mec),      "Smear MC CV with 1.5 #times MEC Model", "l");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_nogtune),  "Smear MC CV with no gTune Model",       "l");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_nopi0tune),"Smear MC CV with #pi^{0} Tune) Model",  "l");
    // leg->AddEntry(h_mcxsec_reco_model.at(k_model_FLUGG), "Smear MC CV with FLUGG  Model",         "l");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_tune1),    "Smear MC CV with Tune 1 Model (Stat.)",         "le");
    leg->Draw();

    // Save and close
    c->Print(Form("plots/run%s/Models/%s/SmearingModelComparison.pdf", _util.run_period, _util.xsec_var));

    fxsec->Close();

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareUnfoldedModels(){

    gStyle->SetOptStat(0);

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH2D* h_temp_2D;
    TH1D* h_temp;

    // The covariance matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_data_cov_tot_unfolded", _util.xsec_var));
    TH2D* h_cov = (TH2D*)h_temp_2D->Clone();
    h_cov->SetDirectory(0);

    // Data XSec
    h_temp  = (TH1D*)fxsec->Get(Form("%s/wiener/h_data_xsec_unfolded", _util.xsec_var));
    TH1D* unf = (TH1D*)h_temp->Clone();
    unf->SetDirectory(0);
    unf->SetLineColor(kBlack);
    unf->SetMarkerStyle(20);
    unf->SetMarkerSize(0.5);
    _util.UndoBinWidthScaling(unf);

    // Additional Smearing matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_ac",_util.xsec_var));
    TH2D* h_ac = (TH2D*)h_temp_2D->Clone();
    h_ac->SetDirectory(0);

    fxsec->Close();


    for (int bin = 1; bin < unf->GetNbinsX()+1; bin++){
        double err = h_cov->GetBinContent(bin, bin);
        unf->SetBinError(bin, std::sqrt(err));
    }
    

    // Now Get the Models
    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Create a vector for the models
    std::vector<std::string> models = {
        "CV",
        "mec",
        "nogtune",
        "nopi0tune",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_CV,
        k_model_mec,
        k_model_nogtune,
        k_model_nopi0tune,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };

    std::vector<TH1D*> h_mcxsec_true_model(models.size());
    std::vector<TH1D*> h_mcxsec_true_model_smear(models.size());

    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
        // MC Xsec True
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", models.at(m).c_str(), vars.at(k_var_trueX).c_str(), vars.at(k_var_trueX).c_str()));
        h_mcxsec_true_model.at(m)         = (TH1D*)h_temp->Clone();
        h_mcxsec_true_model_smear.at(m)   = (TH1D*)h_temp->Clone();

        _util.MatrixMultiply(h_mcxsec_true_model.at(m), h_mcxsec_true_model_smear.at(m), h_ac, "reco_true",false);

    }
    
   
    TLegend *leg = new TLegend(0.4, 0.6, 0.75, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(unf, "Data (Stat. + Sys.)", "ep");
    
    // Now calculate the chi-squared
    double chi, pval;
    int ndof;

    std::cout << "CV" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_CV), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_CV),   Form("MC (CV) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");
    
    std::cout << "1.5 MEC" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_mec), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_mec),   Form("MC (1.5 #times MEC) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");
    
    std::cout << "No gTune" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_nogtune), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_nogtune),   Form("MC (no gTune) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");
    
    std::cout << "no pi0 tune" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_nopi0tune), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_nopi0tune),   Form("MC (no #pi^{0} Tune) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");
    
    // _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_FLUGG), unf, h_cov, chi, ndof, pval);
    // leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_FLUGG),   Form("MC (FLUGG Flux) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");

    std::cout << "Tune 1" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_tune1), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_tune1),   Form("MC (Tune 1) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");

    // Scale the histograms by bin width 
    for (unsigned int m = 0; m < models.size(); m++){
        h_mcxsec_true_model_smear.at(m)->Scale(1.0, "width");
        h_mcxsec_true_model_smear.at(m)->SetLineWidth(2);
    }
    unf->Scale(1.0, "width");

    // Make the plot
    TCanvas *c = new TCanvas("c", "c", 500, 500);
    _util.IncreaseLabelSize( h_mcxsec_true_model_smear.at(k_model_CV), c);
    gPad->SetLeftMargin(0.20);
    c->SetBottomMargin(0.15);
    h_mcxsec_true_model_smear.at(k_model_CV)->GetYaxis()->SetTitleOffset(1.7);
    h_mcxsec_true_model_smear.at(k_model_CV)->SetLineColor(kRed+2);
    
    if (std::string(_util.xsec_var) == "elec_E"){
        h_mcxsec_true_model_smear.at(k_model_CV)->SetMaximum(7);
    }
    else if (std::string(_util.xsec_var) == "elec_ang"){
        h_mcxsec_true_model_smear.at(k_model_CV)->SetMaximum(15);
    }
    else if (std::string(_util.xsec_var) == "elec_cang"){
        h_mcxsec_true_model_smear.at(k_model_CV)->SetMaximum(30.0);
        if (_util.zoom) h_mcxsec_true_model_smear.at(k_model_CV)->GetXaxis()->SetRangeUser(0.6, 1.0);

    }

    h_mcxsec_true_model_smear.at(k_model_CV)->SetMinimum(0.0);
    h_mcxsec_true_model_smear.at(k_model_CV)->Draw("hist");

    h_mcxsec_true_model_smear.at(k_model_mec)->SetLineColor(kGreen+2);
    h_mcxsec_true_model_smear.at(k_model_mec)->Draw("hist,same" );

    h_mcxsec_true_model_smear.at(k_model_nogtune)->SetLineColor(kBlue+2);
    h_mcxsec_true_model_smear.at(k_model_nogtune)->Draw("hist,same" );

    h_mcxsec_true_model_smear.at(k_model_nopi0tune)->SetLineColor(kPink+1);
    h_mcxsec_true_model_smear.at(k_model_nopi0tune)->Draw("hist,same" );

    // h_mcxsec_true_model_smear.at(k_model_FLUGG)->SetLineColor(kYellow+2);
    // h_mcxsec_true_model_smear.at(k_model_FLUGG)->Draw("hist,same" );

    h_mcxsec_true_model_smear.at(k_model_tune1)->SetLineColor(kOrange-1);
    h_mcxsec_true_model_smear.at(k_model_tune1)->Draw("hist,same" );
    
    unf->Draw("E1,X0,same");

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.52, 0.92, 0.52, 0.92);
    

    leg->Draw();
    
    if (_util.zoom)
        c->Print(Form("plots/run%s/Models/%s/DataModelUnfoldedComparison_zoom.pdf", _util.run_period, _util.xsec_var));
    else
        c->Print(Form("plots/run%s/Models/%s/DataModelUnfoldedComparison.pdf", _util.run_period, _util.xsec_var));

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
