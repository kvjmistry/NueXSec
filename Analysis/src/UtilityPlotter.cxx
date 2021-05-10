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
        // StudyPPFXWeights();

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
    // Compare the GENIE and NuWro true Pi0 Energies
    else if (std::string(_util.uplotmode) == "genpi0"){
        // CompareGeneratorPi0();
        CompareSelectedPi0();
        return;
    }
    else if(std::string(_util.uplotmode) == "gencomp"){
        Compare1DFluxGeneratedEvents();
        return;
    }
    // This will call the code to optimise the bin widths
    else if (std::string(_util.uplotmode) == "bins"){
        OptimiseBins();
        return;
    }
    // This will call the function to calculate the covariance matrix
    else if (std::string(_util.uplotmode) == "flux"){
        // CalcFluxCovarianceHP();
        CalcFluxCovarianceBeamline();
        return;
    }
    // This will calculate the event rates for the flux spectrum
    else if (std::string(_util.uplotmode) == "rates"){
        PlotParentEventRates("dk2nu");
        PlotParentEventRates("flugg");
        PlotBeamSimRates();
        return;
    }
    // This will call the function to calculate the covariance matrix
    else if (std::string(_util.uplotmode) == "print"){
        PrintFluxValues();
        PrintXSecResults();
        SaveResponseMatrix();
        return;
    }
    // This will call the code to optimise the bin widths
    else if (std::string(_util.uplotmode) == "models"){ 

        _util.CreateDirectory("Models/" + std::string(_util.xsec_var));
        _util.CreateDirectory("Models/Total");
        TestModelDependence();
        CompareDataCrossSections();
        CompareSmearing();
        CompareUnfoldedModels();
        CompareFakeDataReco();
        CompareFakeDataTrue();
        CompareTotalCrossSec();
        CompareFakeTotalCrossSec();
        CompareTotalDataCrossSections();
        CompareUnfoldedDataCrossSections();
        CheckPi0Coverage();
        CompareMCC8Result();
        ForwardFoldedGeneratorComparison();
        CompareGeneratorTotalCrossSec();
        CompareGeneratorUnfoldedModels();
        CompareXsecPi0Tunings();
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
        nbins = 7;
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

            if (bin == 5 && var == "elec_E"){
                bool fake = true;
                GetFitResult(mean, sigma, 1.43, 3.00, tree, true, fake, true, var);
                lower_bin = 3.00;
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
    _util.CreateDirectory("Resolution/" + std::string(_util.xsec_var));

    _util.CreateDirectory("Purity/" + std::string(_util.xsec_var));

    _util.CreateDirectory("Completeness/" + std::string(_util.xsec_var));


    // Get the vector of bins
    std::vector<double> bins;
    std::string reco_var, true_var;
    
    if (std::string(_util.xsec_var) == "elec_E"){
        bins = _util.reco_shr_bins;
        reco_var = "shr_energy_cali";
        true_var = "elec_e";
    }
    else if (std::string(_util.xsec_var) == "elec_ang"){
        bins = _util.reco_shr_bins_ang;
        reco_var = "effective_angle";
        true_var = "true_effective_angle";
    }
     else if (std::string(_util.xsec_var) == "elec_cang"){
        bins = _util.reco_shr_bins_cang;
        reco_var = "cos_effective_angle";
        true_var = "cos_true_effective_angle";
    }
    
     
    
    // Loop over the bins
    for (float bin = 0; bin < bins.size()-1; bin++ ){

        std::cout <<"\nBin Range: " << bins.at(bin) << " - " << bins.at(bin+1) << " GeV" << std::endl;
        
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "res_reco", reco_var, true_var);
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "res_true", reco_var, true_var);
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "purity", reco_var, true_var);
        PlotQuery(bins.at(bin), bins.at(bin+1), tree, "completeness", reco_var, true_var);


    }
    
}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotQuery(float bin_lower_edge, float bin_upper_edge, TTree* tree, std::string x_var, std::string reco_var, std::string true_var){
    
    TCut query = Form(" weight*(  ((classification.c_str()==\"nue_cc\" || classification.c_str()==\"nuebar_cc\") && passed_selection) && %s > %f && %s < %f )", reco_var.c_str(), bin_lower_edge, reco_var.c_str(), bin_upper_edge); 
    
    TCanvas * c = new TCanvas("c", "c", 500, 500);

    TH1D *htemp;
    if      (x_var == "res_reco" || x_var == "res_true") htemp = new TH1D("htemp","", 30, -1.2, 1.2);
    else if (x_var == "purity") htemp = new TH1D("htemp","", 21, 0, 1.1);
    else if (x_var == "completeness") htemp = new TH1D("htemp","", 21, 0, 1.1);
    else {
        std::cout << "incorrect variable input" << std::endl;
        return;
    }

    
     

    // Draw the Query -- adjust by query type
    if      (x_var == "res_reco") tree->Draw(Form("(%s - %s) / %s >> htemp", reco_var.c_str(), true_var.c_str(), reco_var.c_str()), query);
    else if (x_var == "res_true") tree->Draw(Form("(%s - %s) / %s >> htemp", reco_var.c_str(), true_var.c_str(), true_var.c_str()), query);
    else if (x_var == "purity") tree->Draw("shr_bkt_purity >> htemp", query);
    else if (x_var == "completeness") tree->Draw("shr_bkt_completeness >> htemp", query);
    else {
        std::cout << "incorrect variable input" << std::endl;
        return;
    }    
            
    // Draw the histogram
    htemp->SetLineWidth(2);
    htemp->SetLineColor(kBlack);
    htemp->Draw("hist");

    // Draw the text specifying the bin range
    TLatex* range;
    if (true_var == "elec_e"){
        range  = new TLatex(0.65,0.91, Form("Reco Energy %0.2f - %0.2f GeV",bin_lower_edge, bin_upper_edge ));
    }
    if (true_var == "true_effective_angle"){
        range  = new TLatex(0.65,0.91, Form("Reco #beta %0.2f - %0.2f deg",bin_lower_edge, bin_upper_edge ));
    }
     if (true_var == "cos_true_effective_angle"){
        range  = new TLatex(0.65,0.91, Form("Reco cos#beta %0.2f - %0.2f",bin_lower_edge, bin_upper_edge ));
    }

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

    if (x_var == "res_reco")          htemp->SetTitle("; Reco - True / Reco; Entries");
    else if (x_var == "res_true")     htemp->SetTitle("; Reco - True / True; Entries");
    else if (x_var == "purity")       htemp->SetTitle("; Reco Shower Purity; Entries");
    else if (x_var == "completeness") htemp->SetTitle("; Reco Shower Completeness; Entries");
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
    if (x_var== "res_reco"){
        if (true_var == "elec_e"){
            c->Print(Form("plots/run%s/Resolution/%s/resolution_%0.2fMeV_to_%0.2fMeV_reco.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "true_effective_angle"){
            c->Print(Form("plots/run%s/Resolution/%s/resolution_%0.1fdeg_to_%0.1fdeg_reco.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "cos_true_effective_angle"){
            c->Print(Form("plots/run%s/Resolution/%s/resolution_%0.2f_to_%0.2f_reco.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
    }
    else if (x_var == "res_true"){
        if (true_var == "elec_e"){
            c->Print(Form("plots/run%s/Resolution/%s/resolution_%0.2fMeV_to_%0.2fMeV_true.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "true_effective_angle"){
            c->Print(Form("plots/run%s/Resolution/%s/resolution_%0.1fdeg_to_%0.1fdeg_true.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "cos_true_effective_angle"){
            c->Print(Form("plots/run%s/Resolution/%s/resolution_%0.2f_to_%0.2f_true.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        
    }       
    else if (x_var == "purity"){
        if (true_var == "elec_e"){
            c->Print(Form("plots/run%s/Purity/%s/purity_%0.2fMeV_to_%0.2fMeV.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "true_effective_angle"){
            c->Print(Form("plots/run%s/Purity/%s/purity_%0.1fdeg_to_%0.1fdeg.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "cos_true_effective_angle"){
            c->Print(Form("plots/run%s/Purity/%s/purity_%0.2f_to_%0.2f.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }

        
    }
    else if (x_var == "completeness"){
        if (true_var == "elec_e"){
            c->Print(Form("plots/run%s/Completeness/%s/completeness_%0.2fMeV_to_%0.2fMeV.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "true_effective_angle"){
            c->Print(Form("plots/run%s/Completeness/%s/completeness_%0.1fdeg_to_%0.1fdeg.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        if (true_var == "cos_true_effective_angle"){
            c->Print(Form("plots/run%s/Completeness/%s/completeness_%0.2f_to_%0.2f.pdf", _util.run_period, _util.xsec_var, bin_lower_edge, bin_upper_edge ));
        }
        
    }
    else {
        std::cout << "incorrect variable input" << std::endl;
        return;
    }

    delete htemp;
    delete c;


}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotIntegratedFluxwithThrehold(){
    
    // This changes the plot to average flux or not
    bool draw_averge =true;

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
    h_nue->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu}_{e} / cm^{2} / 5 MeV / 2.0 #times 10^{20} POT");
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
    h_nuebar->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu}_{e} / cm^{2} / 5 MeV / 2.0 #times 10^{20} POT");
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

    summed_flux->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e} + #bar{#nu}_{e} / cm^{2} / 5 MeV / 2.0 #times 10^{20} POT");
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

    // double average_den = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1); // flux denominator nue only
    // double average_den = h_nuebar->Integral(xbin_th, h_nuebar->GetNbinsX()+1); // flux denominator nuebar only

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
    h_summed_flux_clone->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e} + #bar{#nu}_{e} / cm^{2} / 2.0 #times 10^{20} POT");
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
    _util.Draw_ubooneSim(c, 0.37, 0.93, 0.37, 0.91);

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
    _util.GetFile(f_mc, "../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root"); // Get the run 1 MC file
    _util.GetTree(f_mc, mc_tree, "nuselection/NeutrinoSelectionFilter");

    SliceContainer SC;
    SC.Initialise(mc_tree, _util.k_mc, _util);

    std::vector<std::string> vars = {"nu_e", "elec_e"};

    // Create histograms for the hit purity
    std::vector<std::vector<TH2D*>> h_hit_pur;
    
    // 1D pi0 momentum
    TH1D *h_pi0_momentum = new TH1D("h_true_pi0_momentum", "; #pi^{0} Momentum [GeV/c]; Entries", 40, 0, 2.0);
    
    // 2D shower multiplicity vd nue/electron energy
    TH2D *h_shr_multi_nue_E         = new TH2D("h_shr_multi_nue_E", "; Shower Multiplicty;#nu_{e} Energy [GeV] ", 6, 0, 6, 10, 0, 4.0);
    TH2D *h_shr_multi_elec_e        = new TH2D("h_shr_multi_elec_e", "; Shower Multiplicty;Electron Energy [GeV] ", 6, 0, 6, 10, 0, 4.0);
    TH2D *h_shr_multi_nuebar_E      = new TH2D("h_shr_multi_nuebar_E", "; Shower Multiplicty;#bar{#nu}_{e} Energy [GeV] ", 6, 0, 6, 10, 0, 4.0);
    TH2D *h_shr_multi_elec_e_nuebar = new TH2D("h_shr_multi_elec_e_nuebar", "; Shower Multiplicty;Positron Energy [GeV] ", 6, 0, 6, 10, 0, 4.0);
    
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

        // Classify the event -- sets variable in the slice contianer
        SC.SliceClassifier(_util.k_mc);      // Classification of the event

        // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
        SC.SetThresholdEvent(_util.k_mc);

        // If the backtracked pdg of the leading shower is not an electron then alter classification
        SC.SetNonLdgShrEvent(_util.k_mc);
        
        SC.SliceInteractionType(_util.k_mc); // Genie interaction type
        SC.ParticleClassifier(_util.k_mc);   // The truth matched particle type of the leading shower
        SC.Pi0Classifier(_util.k_mc); 

        // Set derived variables in the slice container
        SC.SetSignal();                // Set the event as either signal or other
        SC.SetFakeData();              // Set the classifcation as data if fake data mode
        SC.SetTrueElectronThetaPhi();  // Set the true electron theta and phi variables
        SC.SetNuMIAngularVariables();  // Set the NuMI angular variables
        SC.CalibrateShowerEnergy();    // Divide the shower energy by 0.83 so it is done in one place

        bool is_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 

        double weight = _util.GetCVWeight(_util.k_mc, SC.weightSplineTimesTune, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv, SC.interaction, SC.elec_e);
        
        // True nue energy
        h_hit_pur.at(0).at(SC.classification.second)->Fill(SC.nu_e, SC.nu_purity_from_pfp, weight);
        
        // True electron energy
        h_hit_pur.at(1).at(SC.classification.second)->Fill(SC.elec_e, SC.nu_purity_from_pfp, weight);

        // Pi0 Momentum
        if (is_in_fv) 
            h_pi0_momentum->Fill(std::sqrt(SC.pi0_e*SC.pi0_e - 0.134*0.134), weight);

        // Nue cc
        if (SC.nslice == 1 && SC.nu_pdg == 12 && is_in_fv && SC.nu_purity_from_pfp > 0.5 && SC.ccnc == _util.k_CC){

            h_shr_multi_nue_E->Fill(SC.n_showers, SC.nu_e, weight);
            h_shr_multi_elec_e->Fill(SC.n_showers, SC.elec_e, weight);
        }
        // nuebar cc
        if (SC.nslice == 1 && SC.nu_pdg == -12 && is_in_fv && SC.nu_purity_from_pfp > 0.5 && SC.ccnc == _util.k_CC){
            h_shr_multi_nuebar_E->Fill(SC.n_showers, SC.nu_e, weight);
            h_shr_multi_elec_e_nuebar->Fill(SC.n_showers, SC.elec_e, weight);
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

    RowNorm(h_shr_multi_nue_E);
    RowNorm(h_shr_multi_elec_e);
    RowNorm(h_shr_multi_nuebar_E);
    RowNorm(h_shr_multi_elec_e_nuebar);

    // h_shr_multi_nue_E->GetZaxis()->SetRangeUser(0, 0.35);
    // h_shr_multi_nuebar_E->GetZaxis()->SetRangeUser(0, 0.35);
    // h_shr_multi_elec_e_nuebar->GetZaxis()->SetRangeUser(0, 0.65);
    // h_shr_multi_elec_e->GetZaxis()->SetRangeUser(0, 0.65);

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
void UtilityPlotter::RowNorm(TH2D* hist){

    // Loop over rows
    for (int col=1; col<hist->GetYaxis()->GetNbins()+1; col++) {
        double integral = 0;

        // Loop over columns and get the integral
        for (int row=1; row<hist->GetXaxis()->GetNbins()+1; row++){
            integral+=hist->GetBinContent(row, col);            
        }

        // Now normalise the column entries by the integral
        for (int row=1; row<hist->GetXaxis()->GetNbins()+1; row++){
            hist->SetBinContent(row, col, hist->GetBinContent(row, col)/ integral );
            
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
        h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run1_%s_0_smearing", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), models.at(m).c_str()));
        if (h_temp_2D == NULL) std::cout <<"Help!" << m << std::endl;
        h_response_model.at(m) = (TH2D*)h_temp_2D->Clone();

        // MC xsec in True
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), _util.vars.at(k_var_trueX).c_str()));
        h_mcxsec_true_model.at(m) = (TH1D*)h_temp->Clone();
        h_mcxsec_reco_model.at(m) = (TH1D*)h_dataxsec->Clone();
       
        _util.MatrixMultiply(h_mcxsec_true_model.at(m), h_mcxsec_reco_model.at(m), h_response_model.at(k_model_CV), "true_reco", true);
    }

    // Set the line colours
    h_mcxsec_reco_model.at(k_model_CV)       ->SetLineColor(kRed+2);
    h_mcxsec_reco_model.at(k_model_mec)      ->SetLineColor(kGreen+2);
    h_mcxsec_reco_model.at(k_model_nogtune)  ->SetLineColor(kBlue+2);
    h_mcxsec_reco_model.at(k_model_nopi0tune)->SetLineColor(kPink+1);
    h_mcxsec_reco_model.at(k_model_FLUGG)    ->SetLineColor(kViolet-1);
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
    h_mcxsec_reco_model.at(k_model_FLUGG)->Draw("hist,same");
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
    
    
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_FLUGG), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_FLUGG),  Form("MC (FLUGG) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");

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
        c->Print(Form("plots/run%s/Models/%s/run%s_DataModelComparison_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));

    fxsec->Close();

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareSmearing(){

    gStyle->SetOptStat(0);

    // Load in the cross section output
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
        "geniev3",
        "nuwro",
        "FLUGG",
        // "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_mec,
        k_model_geniev3,
        k_model_nuwro,
        k_model_FLUGG,
        // k_model_tune1,
        k_MODEL_MAX
    };

    // Create the vector of histograms
    std::vector<TH2D*> h_response_model(models.size());
    std::vector<TH1D*> h_mcxsec_reco_model(models.size());
    
    // MC Xsec True
    h_temp  = (TH1D*)fxsec->Get(Form("CV/%s/h_run1_CV_0_%s_mc_xsec",_util.vars.at(k_var_trueX).c_str(), _util.vars.at(k_var_trueX).c_str()));
    TH1D* h_mcxsec_true         = (TH1D*)h_temp->Clone();
    
    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
        
        // clone to get the binning, these hists get reset
        h_mcxsec_reco_model.at(m) = (TH1D*)h_mcxsec_reco->Clone();

        // Get the response matrix
        h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run1_%s_0_smearing", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), models.at(m).c_str()));
        h_response_model.at(m) = (TH2D*)h_temp_2D->Clone();

        // Apply the response matrix to the CV MC True dist
        _util.MatrixMultiply(h_mcxsec_true, h_mcxsec_reco_model.at(m), h_response_model.at(m), "true_reco", true);
    }

    // Set the line colours
    h_mcxsec_reco_model.at(k_model_mec)      ->SetLineColor(kGreen+2);
    h_mcxsec_reco_model.at(k_model_geniev3)  ->SetLineColor(kBlue+2);
    h_mcxsec_reco_model.at(k_model_nuwro)->SetLineColor(kPink+1);
    h_mcxsec_reco_model.at(k_model_FLUGG)    ->SetLineColor(kViolet-1);
    // h_mcxsec_reco_model.at(k_model_tune1)    ->SetLineColor(kOrange-1);

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
    h_mcxsec_reco_model.at(k_model_geniev3)  ->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_nuwro)    ->Draw("hist,E,same");
    h_mcxsec_reco_model.at(k_model_FLUGG)    ->Draw("hist,E,same");
    // h_mcxsec_reco_model.at(k_model_tune1)    ->Draw("hist,E,same");
    h_mcxsec_reco->Draw("hist,E,same");

    // Create the legend
    TLegend *leg = new TLegend(0.4, 0.5, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_mcxsec_reco, "MC (Stat.)", "el");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_mec),      "Smear MC CV with 1.5 #times MEC Model", "l");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_geniev3),  "Smear MC CV with GENIE v3 Model",       "l");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_nuwro),    "Smear MC CV with NuWro Model",   "le");
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_FLUGG),    "Smear MC CV with FLUGG  Model (Stat.)",         "le");
    // leg->AddEntry(h_mcxsec_reco_model.at(k_model_tune1),    "Smear MC CV with Tune 1 Model (Stat.)", "le");
    leg->Draw();

    // Save and close
    c->Print(Form("plots/run%s/Models/%s/run%s_SmearingModelComparison_%s.pdf", _util.run_period, _util.xsec_var,  _util.run_period, _util.xsec_var));

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
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), _util.vars.at(k_var_trueX).c_str()));
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
    
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_FLUGG), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_FLUGG),   Form("MC (FLUGG) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");

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

    h_mcxsec_true_model_smear.at(k_model_FLUGG)->SetLineColor(kViolet-1);
    h_mcxsec_true_model_smear.at(k_model_FLUGG)->Draw("hist,same" );

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
        c->Print(Form("plots/run%s/Models/%s/run%s_DataModelUnfoldedComparison_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareFakeDataReco(){

    gStyle->SetOptStat(0);

    /// Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH2D* h_temp_2D;

    // Get the MC covariance Matrix
    h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/er/h_cov_sys_mcxsec_reco", _util.xsec_var ));
    TH2D* h_cov_tot = (TH2D*)h_temp_2D->Clone();
    h_cov_tot->SetDirectory(0);

    h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/er/h_cov_xsec_sys_mcxsec_reco", _util.xsec_var ));
    TH2D* h_cov_xsec = (TH2D*)h_temp_2D->Clone();
    h_cov_xsec->SetDirectory(0);

    h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/er/h_cov_flux_sys_mcxsec_reco", _util.xsec_var ));
    TH2D* h_cov_flux = (TH2D*)h_temp_2D->Clone();
    h_cov_flux->SetDirectory(0);

    fxsec->Close();

    /// Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Create a vector for the models
    std::vector<std::string> models = {
        "Input",
        "mec",
        "geniev3",
        "nuwro",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_input,
        k_model_mec,
        k_model_geniev3,
        k_model_nuwro,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };
    
    // Temp histograms
    TH1D* htemp;
    TH2D* htemp2D;


    std::vector<TH1D*> h_true(k_MODEL_MAX);
    std::vector<TH1D*> h_true_smear(k_MODEL_MAX);
    std::vector<TH1D*> h_fake(k_MODEL_MAX);
    std::vector<TH2D*> h_response(k_MODEL_MAX);
    std::vector<TH2D*> h_cov_m(k_MODEL_MAX);
    
    double ymax = 1.0;

    // Get the cv hist 
    htemp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", "CV", _util.vars.at(k_var_recoX).c_str(), _util.vars.at(k_var_recoX).c_str()));
    TH1D *h_temp_CV = (TH1D*)htemp->Clone();
    h_temp_CV->Scale(1.0, "width");
    
    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){

        htemp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec",models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), _util.vars.at(k_var_trueX).c_str()));
        h_true.at(m)        = (TH1D*)htemp->Clone();
        h_true_smear.at(m)  = (TH1D*)htemp->Clone();

        htemp  = (TH1D*)fxsec->Get(Form("fake%s/%s/h_run1_CV_0_%s_data_xsec",models.at(m).c_str(), _util.vars.at(k_var_recoX).c_str(), _util.vars.at(k_var_recoX).c_str()));
        h_fake.at(m)        = (TH1D*)htemp->Clone();

        // Get the response matrix
        htemp2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run1_%s_0_smearing", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), models.at(m).c_str()));
        h_response.at(m) = (TH2D*)htemp2D->Clone();

        // Get the Covariance matrix
        if (m == k_model_FLUGG)
            h_cov_m.at(m) = (TH2D*) h_cov_flux->Clone();
        else
            h_cov_m.at(m) = (TH2D*) h_cov_xsec->Clone();

        // Apply the response matrix to the model MC True dist to get the reco dist back
        _util.MatrixMultiply(h_true.at(m), h_true_smear.at(m), h_response.at(m), "true_reco", true);

        // Set the line colours
        h_true_smear.at(m)->SetLineColor(kRed+2);
        h_fake.at(m)->SetLineColor(kBlack);
        h_true_smear.at(m)->SetLineWidth(2);
        // h_fake.at(m)->SetLineWidth(2);
        h_fake.at(m)->SetMarkerStyle(20);
        h_fake.at(m)->SetMarkerSize(0.5);

        // h_true_smear.at(m)->Scale(1.0, "width");
        h_fake.at(m)->Scale(1.0, "width");

        // Convert the Covariance Matrix-- switching from MC CV deviations to Fake Data CV deviation
        _util.ConvertCovarianceUnits(h_cov_m.at(m),
                               h_temp_CV,
                               h_fake.at(m));

        // Now set the bin errors
        for (int bin = 1; bin <  h_fake.at(m)->GetNbinsX()+1; bin++ ){    
            h_fake.at(m)->SetBinError(bin, std::sqrt(h_cov_m.at(m)->GetBinContent(bin, bin)));
        }

        ymax = h_true_smear.at(m)->GetMaximum();

        if (h_fake.at(m)->GetMaximum() > h_true_smear.at(m)->GetMaximum())
            ymax = h_fake.at(m)->GetMaximum();
    }

    // Now lets plot
    TCanvas *c;
    
    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
    
        c = new TCanvas("c", "c", 500, 500);
        c->SetLeftMargin(0.2);
        c->SetBottomMargin(0.15);
        h_true_smear.at(m)->GetYaxis()->SetTitleOffset(1.7);
        h_true_smear.at(m)->SetMinimum(0);
        h_true_smear.at(m)->SetMaximum(ymax + 0.4*ymax);

        // Set the line colours
        if (m == k_model_input)    h_true_smear.at(k_model_input)    ->SetLineColor(kRed+2);
        if (m == k_model_mec)      h_true_smear.at(k_model_mec)      ->SetLineColor(kGreen+2);
        if (m == k_model_geniev3)  h_true_smear.at(k_model_geniev3)  ->SetLineColor(kBlue+2);
        if (m == k_model_nuwro)    h_true_smear.at(k_model_nuwro)    ->SetLineColor(kPink+1);
        if (m == k_model_FLUGG)    h_true_smear.at(k_model_FLUGG)    ->SetLineColor(kViolet-1);
        if (m == k_model_tune1)    h_true_smear.at(k_model_tune1)    ->SetLineColor(kOrange-1);

        h_true_smear.at(m)->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());

        TH1D* h_error_hist = (TH1D*)h_true_smear.at(m)->Clone();
        h_error_hist->SetFillColorAlpha(12, 0.15);

        h_true_smear.at(m)->Draw("hist");

        if (m != k_model_input)
            h_true_smear.at(k_model_input)->Draw("hist,same");

        h_error_hist->Draw("e2, same");
        h_fake.at(m)->Draw("E,same");

        // Create the legend
        TLegend *leg = new TLegend(0.4, 0.5, 0.85, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        if (m == k_model_input) models.at(m) = "CV";

        if (m != k_model_input)
            leg->AddEntry( h_true_smear.at(k_model_input), "CV", "lf");

        leg->AddEntry(h_error_hist, Form("True %s (stat.)", models.at(m).c_str()), "lf");
        leg->AddEntry(h_fake.at(m), Form("Fake %s (sys.)", models.at(m).c_str()), "elp");
        leg->Draw();

        // Save and close
        c->Print(Form("plots/run%s/Models/%s/run%s_FakeDataComparison_%s_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var, models.at(m).c_str()));
        delete c;
    }

    fxsec->Close();

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareFakeDataTrue(){

    gStyle->SetOptStat(0);

    // Create a vector for the models
    std::vector<std::string> models = {
        "Input",
        "mec",
        "geniev3",
        "nuwro",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_input,
        k_model_mec,
        k_model_geniev3,
        k_model_nuwro,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };


    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH2D* h_temp_2D;
    TH1D* h_temp;

    // True MC Xsec
    h_temp  = (TH1D*)fxsec->Get(Form("%s/wiener/h_mc_xsec_true", _util.xsec_var));
    TH1D* h_true_mc_xsec = (TH1D*)h_temp->Clone();
    h_true_mc_xsec->SetDirectory(0);

    // Response Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_response",_util.xsec_var));
    TH2D* h_response = (TH2D*)h_temp_2D->Clone();
    h_response->SetDirectory(0);

    // Total Covariance Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_cov_sys_dataxsec_reco",_util.xsec_var));
    TH2D* h_cov_reco = (TH2D*)h_temp_2D->Clone();
    h_cov_reco->SetDirectory(0);

    // Flux Covariance Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_cov_flux_sys_dataxsec_reco",_util.xsec_var));
    TH2D* h_cov_reco_flux = (TH2D*)h_temp_2D->Clone();
    h_cov_reco_flux->SetDirectory(0);

    // Xsec Covariance Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_cov_xsec_sys_dataxsec_reco",_util.xsec_var));
    TH2D* h_cov_reco_xsec = (TH2D*)h_temp_2D->Clone();
    h_cov_reco_xsec->SetDirectory(0);

    fxsec->Close();

    // Now Get the Models
    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Load in the CV data cross section
    h_temp  = (TH1D*)fxsec->Get(Form("CV/%s/h_run%s_CV_0_%s_data_xsec",_util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
    TH1D* h_reco_data_xsec = (TH1D*)h_temp->Clone();

    std::vector<TH1D*> h_true(k_MODEL_MAX);
    std::vector<TH1D*> h_fake(k_MODEL_MAX);
    std::vector<TH2D*> h_cov_diag(k_MODEL_MAX);

    TH1D* h_fake_xsec_smear_CV;

    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){

        // Get true model xsec
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run%s_CV_0_%s_mc_xsec", models.at(m).c_str(),_util.vars.at(k_var_trueX).c_str(), _util.run_period, _util.vars.at(k_var_trueX).c_str()));
        h_true.at(m) = (TH1D*)h_temp->Clone();

        // Get fake Tune1 data
        h_temp  = (TH1D*)fxsec->Get(Form("fake%s/%s/h_run%s_CV_0_%s_data_xsec", models.at(m).c_str(), _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
        h_fake.at(m) = (TH1D*)h_temp->Clone();


        // Set diagonals of covariance matrix
        if (m == k_model_FLUGG)
            h_cov_diag.at(m) = (TH2D*)h_cov_reco_flux->Clone();
        else
            h_cov_diag.at(m) = (TH2D*)h_cov_reco_xsec->Clone();
        
        // Convert the Covariance Matrix-- switching from data deviations to Fake Data CV deviation
        _util.ConvertCovarianceUnits(h_cov_diag.at(m),
                               h_reco_data_xsec,
                               h_fake.at(m));


        // Initialise the WienerSVD class
        WienerSVD _wSVD;
        _wSVD.Initialise(_util);
        _wSVD.DoUnfolding(2, 0, h_true_mc_xsec, h_fake.at(m), h_response, h_cov_diag.at(m));

        for (int bin = 1; bin < _wSVD.unf->GetNbinsX()+1; bin++){
            double err = _wSVD.unfcov->GetBinContent(bin, bin);
            _wSVD.unf->SetBinError(bin, std::sqrt(err));
        }
        
        _wSVD.unf->SetLineColor(kBlack);
        _wSVD.unf->SetMarkerStyle(20);
        _wSVD.unf->SetMarkerSize(0.5);


        // Matrix multiply the true xsec by AC
        TH1D* h_fake_xsec_smear = (TH1D*)h_true.at(m)->Clone();
        _util.MatrixMultiply(h_true.at(m), h_fake_xsec_smear, _wSVD.smear, "reco_true", false);
        
        // Make the plot
        TCanvas *c = new TCanvas("c", "c", 500, 500);
        // _util.IncreaseLabelSize( h_mcxsec_true_model_smear.at(k_model_CV), c);
        gPad->SetLeftMargin(0.20);
        c->SetBottomMargin(0.15);


        h_fake_xsec_smear->Scale(1.0, "width");
        _wSVD.unf->Scale(1.0, "width");

        double ymax = h_fake_xsec_smear->GetMaximum();

        if (_wSVD.unf->GetMaximum() > ymax)
            ymax = _wSVD.unf->GetMaximum();

        h_fake_xsec_smear->SetMaximum(ymax + ymax*0.4);
        h_fake_xsec_smear->SetLineWidth(2);
        h_fake_xsec_smear->SetLineColor(kRed+2);

        if (std::string(_util.xsec_var) == "elec_E"){
            h_fake_xsec_smear->SetMaximum(8);
        }
        else if (std::string(_util.xsec_var) == "elec_ang"){
            h_fake_xsec_smear->SetMaximum(15);
        }
        else if (std::string(_util.xsec_var) == "elec_cang"){
            h_fake_xsec_smear->SetMaximum(30.0);
        }

        h_fake_xsec_smear->SetMinimum(0);

        // Set the line colours
        if (m == k_model_input)    h_fake_xsec_smear->SetLineColor(kRed+2);
        if (m == k_model_mec)      h_fake_xsec_smear->SetLineColor(kGreen+2);
        if (m == k_model_geniev3)  h_fake_xsec_smear->SetLineColor(kBlue+2);
        if (m == k_model_nuwro)    h_fake_xsec_smear->SetLineColor(kPink+1);
        if (m == k_model_FLUGG)    h_fake_xsec_smear->SetLineColor(kViolet-1);
        if (m == k_model_tune1)    h_fake_xsec_smear->SetLineColor(kOrange-1);

        if (m == k_model_input)
            h_fake_xsec_smear_CV = (TH1D*)h_fake_xsec_smear->Clone();

        TH1D* h_error_hist = (TH1D*)h_fake_xsec_smear->Clone();
        h_error_hist->SetFillColorAlpha(12, 0.15);

        h_fake_xsec_smear->Draw("hist");

        if (m != k_model_input)
            h_fake_xsec_smear_CV->Draw("hist,same");

        h_error_hist->Draw("E2, same");
        _wSVD.unf->Draw("E,same");


        // Create the legend
        TLegend *leg = new TLegend(0.4, 0.5, 0.85, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        if (m == k_model_input) models.at(m) = "CV";

        if (m != k_model_input)
            leg->AddEntry(h_fake_xsec_smear_CV, "CV", "lf");
        leg->AddEntry(h_error_hist, Form("True %s (stat.)", models.at(m).c_str()), "lf");
        leg->AddEntry(_wSVD.unf, Form("Fake %s (sys.)", models.at(m).c_str()), "elp");
        leg->Draw();

        
        c->Print(Form("plots/run%s/Models/%s/run%s_UnfoldedFakeDataComparison_%s_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var, models.at(m).c_str()));
        delete c;

        delete _wSVD.smear;
        delete _wSVD.wiener;
        delete _wSVD.unfcov;
        delete _wSVD.unf;
        delete _wSVD.diff;
        delete _wSVD.bias;
        delete _wSVD.bias2;
        delete _wSVD.fracError;
        delete _wSVD.absError;
        delete _wSVD.MSE;
        delete _wSVD.MSE2;
    }

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareTotalCrossSec(){

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c", 150, 350);
    gPad->SetLeftMargin(0.3);

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH1D* h_temp;

    // Data Xsec stat
    h_temp  = (TH1D*)fxsec->Get("total/h_data_xsec_stat_reco");
    TH1D* h_data_stat = (TH1D*)h_temp->Clone();
    h_data_stat->SetDirectory(0);

    // Data Xsec stat + sys
    h_temp  = (TH1D*)fxsec->Get("total/h_data_xsec_stat_sys_reco");
    TH1D* h_data = (TH1D*)h_temp->Clone();
    h_data->SetDirectory(0);

    fxsec->Close();

    // Now Get the Models
    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Now get some other models
    // Create a vector for the models
    std::vector<std::string> models = {
        "CV",
        "mec",
        "nogtune",
        "nuwro",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_CV,
        k_model_mec,
        k_model_nogtune,
        k_model_nuwro,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };

    std::vector<TH1D*> h_model_xsec(k_MODEL_MAX);
    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){

        // Get true tune1 xsec
        h_temp  = (TH1D*)fxsec->Get(Form("%s/integrated/h_run%s_CV_0_integrated_mc_xsec", models.at(m).c_str(), _util.run_period));
        h_model_xsec.at(m) = (TH1D*)h_temp->Clone();
        h_model_xsec.at(m)->SetLineWidth(2);
       
    }


    h_model_xsec.at(k_model_CV)->SetLineColor(kRed+2);
    h_model_xsec.at(k_model_mec)->SetLineColor(kGreen+2);
    h_model_xsec.at(k_model_nogtune)->SetLineColor(kBlue+2);
    h_model_xsec.at(k_model_nuwro)->SetLineColor(kPink+1);
    h_model_xsec.at(k_model_FLUGG)->SetLineColor(kViolet-1);
    h_model_xsec.at(k_model_tune1)->SetLineColor(kOrange-1);


    // X-Axis
    h_data->GetXaxis()->SetRangeUser(0.0,1.0); 
    h_data->GetXaxis()->SetLabelOffset(999);
    h_data->GetXaxis()->SetLabelSize(0);
    h_data->GetXaxis()->SetTickLength(0);
    
    // Y-Axis
    h_data->GetYaxis()->SetRangeUser(3.0, 8.0);
    h_data->GetYaxis()->CenterTitle();
    h_data->GetYaxis()->SetLabelSize(0.1);
    h_data->GetYaxis()->SetTitleSize(0.1);
   
    // h_data->SetLineWidth(2);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.7);
    h_data->SetLineColor(kBlack);
    h_data->Draw("E1,X0");

    // Statistical band
    h_data_stat->SetLineColor(kBlack);
    h_data_stat->Draw("E1,X0,same");

    // Draw the models
    for (unsigned int m = 0; m < models.size(); m++){
        h_model_xsec.at(m)->Draw("hist,same");
    }


    h_data->Draw("E1,X0,same");
    h_data_stat->Draw("E1,X0,same");


    // Draw the Legend
    TLegend *leg = new TLegend(0.35, 0.70, 0.70, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    h_data->SetMarkerSize(0.4);
    leg->AddEntry(h_data, "Data (stat. + sys.)",        "ep");
    leg->AddEntry( h_model_xsec.at(k_model_CV),         "MC (CV)", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_mec),        "MC (1.5 #times MEC)", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_nogtune),    "MC (no gTune)", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_nuwro),      "MC (NuWro)", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_FLUGG),      "MC (FLUGG Flux)", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_tune1),      "MC (Tune 1)", "lf");
    
    leg->Draw();

    gStyle->SetLegendTextSize(0.06);

    double Data_POT; 

    // Set the scale factors
    if (strcmp(_util.run_period, "1") == 0){
        Data_POT = _util.config_v.at(_util.k_Run1_Data_POT); // Define this variable here for easier reading
    }
    else if (strcmp(_util.run_period, "3") == 0){
        Data_POT = _util.config_v.at(_util.k_Run3_Data_POT); // Define this variable here for easier reading
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
        exit(1);
    }

    Data_POT = Data_POT / 1.0e20;

    TLatex *t = new TLatex(.34, .145, Form("#splitline{MicroBooNE NuMI}{Data %2.1f#times10^{20} POT}", Data_POT));
    t->SetTextColor(kBlack);
    t->SetNDC();
    t->SetTextSize(2.0/30.);
    t->SetTextAlign(11);
    t->Draw();


    c->Print(Form("plots/run%s/Models/Total/run%s_TotalCrossSectionComparison.pdf", _util.run_period, _util.run_period));

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareFakeTotalCrossSec(){

    gStyle->SetOptStat(0);

    // Now Get the Models
    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Now get some other models
    // Create a vector for the models
    std::vector<std::string> models = {
        "Input",
        "mec",
        "geniev3",
        "nuwro",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_input,
        k_model_mec,
        k_model_geniev3,
        k_model_nuwro,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };

    TH1D* h_temp;

    std::vector<TH1D*> h_model_xsec(k_MODEL_MAX);
    std::vector<TH1D*> h_gen(k_MODEL_MAX);
    std::vector<TH1D*> h_fake_xsec(k_MODEL_MAX);
    
    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){

        TCanvas *c = new TCanvas("c", "c", 150, 350);
        gPad->SetLeftMargin(0.3);

        // Get true tune1 xsec
        h_temp  = (TH1D*)fxsec->Get(Form("%s/integrated/h_run%s_CV_0_integrated_mc_xsec", models.at(m).c_str(), _util.run_period));
        h_model_xsec.at(m) = (TH1D*)h_temp->Clone();
        h_model_xsec.at(m)->SetLineWidth(2);
        h_model_xsec.at(m)->SetLineColor(kRed+2);

        h_temp  = (TH1D*)fxsec->Get(Form("%s/integrated/h_run%s_CV_0_integrated_gen", models.at(m).c_str(), _util.run_period));
        h_gen.at(m) = (TH1D*)h_temp->Clone();

        // Set the bin error to use the generated events err
        h_model_xsec.at(m)->SetBinError(1,  h_model_xsec.at(m)->GetBinContent(1) * h_gen.at(m)->GetBinError(1) / h_gen.at(m)->GetBinContent(1));


        // Get total mc xsec fake data prediction
        h_temp  = (TH1D*)fxsec->Get(Form("fake%s/integrated/h_run%s_CV_0_integrated_data_xsec", models.at(m).c_str(), _util.run_period));
        h_fake_xsec.at(m) = (TH1D*)h_temp->Clone();

        // Set the error to be equal to the total systematic uncertainty of ~21%
        if (m == k_model_FLUGG)
            h_fake_xsec.at(m)->SetBinError(1,h_fake_xsec.at(m)->GetBinContent(1) * 0.216 );
        else
            h_fake_xsec.at(m)->SetBinError(1,h_fake_xsec.at(m)->GetBinContent(1) * std::sqrt(0.032*0.032 + 0.056*0.056) ); 
            // h_fake_xsec.at(m)->SetBinError(1,h_fake_xsec.at(m)->GetBinContent(1) * std::sqrt(0.04551*0.04551 + 0.0374*0.0374) );

        // Set the line colours
        if (m == k_model_input)    h_model_xsec.at(k_model_input)    ->SetLineColor(kRed+2);
        if (m == k_model_mec)      h_model_xsec.at(k_model_mec)      ->SetLineColor(kGreen+2);
        if (m == k_model_geniev3)  h_model_xsec.at(k_model_geniev3)  ->SetLineColor(kBlue+2);
        if (m == k_model_nuwro)    h_model_xsec.at(k_model_nuwro)    ->SetLineColor(kPink+1);
        if (m == k_model_FLUGG)    h_model_xsec.at(k_model_FLUGG)    ->SetLineColor(kViolet-1);
        if (m == k_model_tune1)    h_model_xsec.at(k_model_tune1)    ->SetLineColor(kOrange-1);
        

        // X-Axis
        h_fake_xsec.at(m)->GetXaxis()->SetRangeUser(0.0,1.0); 
        h_fake_xsec.at(m)->GetXaxis()->SetLabelOffset(999);
        h_fake_xsec.at(m)->GetXaxis()->SetLabelSize(0);
        h_fake_xsec.at(m)->GetXaxis()->SetTickLength(0);
        
        // Y-Axis
        h_fake_xsec.at(m)->GetYaxis()->SetRangeUser(3.0, 9.0);
        h_fake_xsec.at(m)->GetYaxis()->CenterTitle();
        h_fake_xsec.at(m)->GetYaxis()->SetLabelSize(0.1);
        h_fake_xsec.at(m)->GetYaxis()->SetTitleSize(0.1);

        h_fake_xsec.at(m)->SetMarkerStyle(20);
        h_fake_xsec.at(m)->SetMarkerSize(0.7);
        h_fake_xsec.at(m)->SetLineColor(kBlack);
        h_fake_xsec.at(m)->Draw("E1,X0");
        h_model_xsec.at(m)->Draw("hist,E,same");
        h_fake_xsec.at(m)->Draw("E1,X0,same");

        if (m != k_model_input)
            h_model_xsec.at(k_model_input)->Draw("hist,E,same");


        // Draw the Legend
        TLegend *leg = new TLegend(0.35, 0.70, 0.70, 0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        h_fake_xsec.at(m)->SetMarkerSize(0.4);
        if (m == k_model_input) models.at(m) = "CV";

        if (m != k_model_input)
            leg->AddEntry(h_model_xsec.at(k_model_input), "CV", "l");

        leg->AddEntry(h_model_xsec.at(m), Form("True %s", models.at(m).c_str()), "l");
        leg->AddEntry(h_fake_xsec.at(m),  Form("Fake %s", models.at(m).c_str()),  "ep");
        leg->Draw();

        gStyle->SetLegendTextSize(0.06);

        c->Print(Form("plots/run%s/Models/Total/run%s_FakeTotalCrossSectionComparison_%s.pdf", _util.run_period, _util.run_period,  models.at(m).c_str()));
        delete c;
    
    }
}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareDataCrossSections(){

    gStyle->SetOptStat(0);

    // std::string genmode = "flux";
    std::string genmode = "other";

    /// Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH1D* h_temp;

    // Data XSec
    h_temp  = (TH1D*)fxsec->Get(Form("%s/er/h_data_xsec_sys_reco", _util.xsec_var));
    TH1D* h_dataxsec = (TH1D*)h_temp->Clone();
    h_dataxsec->SetDirectory(0);
    h_dataxsec->SetLineColor(kBlack);
    h_dataxsec->SetMarkerStyle(20);
    h_dataxsec->SetMarkerSize(0.5);

    // Get the MC covariance Matrix
    TH2D* h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/er/h_cov_sys_mcxsec_reco", _util.xsec_var ));
    TH2D* h_cov_tot = (TH2D*)h_temp_2D->Clone();
    h_cov_tot->SetDirectory(0);

    h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/er/h_cov_xsec_sys_mcxsec_reco", _util.xsec_var ));
    TH2D* h_cov_xsec = (TH2D*)h_temp_2D->Clone();
    h_cov_xsec->SetDirectory(0);

    h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/er/h_cov_flux_sys_mcxsec_reco", _util.xsec_var ));
    TH2D* h_cov_flux = (TH2D*)h_temp_2D->Clone();
    h_cov_flux->SetDirectory(0);

    fxsec->Close();

    TH2D* h_cov_reco;

    // Create a vector for the models
    std::vector<std::string> models = {
        "Input",
        "mec",
        "geniev3",
        "nuwro",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_input,
        k_model_mec,
        k_model_geniev3,
        k_model_nuwro,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };

    std::vector<TH1D*> h_dataxsec_model(models.size());

    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
    
        // Data X Sec MEC
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_data_xsec", models.at(m).c_str(), _util.vars.at(k_var_recoX).c_str(), _util.vars.at(k_var_recoX).c_str()));
        h_dataxsec_model.at(m) = (TH1D*)h_temp->Clone();
        h_dataxsec_model.at(m)->Scale(1.0, "width");
        h_dataxsec_model.at(m)->SetLineWidth(2);

        // Set the line colours
        if (m == k_model_mec)      h_dataxsec_model.at(k_model_mec)      ->SetLineColor(kGreen+2);
        if (m == k_model_geniev3)  h_dataxsec_model.at(k_model_geniev3)  ->SetLineColor(kBlue+2);
        if (m == k_model_nuwro)    h_dataxsec_model.at(k_model_nuwro)    ->SetLineColor(kPink+1);
        if (m == k_model_FLUGG)    h_dataxsec_model.at(k_model_FLUGG)    ->SetLineColor(kViolet-1);
        if (m == k_model_tune1)    h_dataxsec_model.at(k_model_tune1)    ->SetLineColor(kOrange-1);
    }

    // Get the Covariance matrix
    if (genmode == "flux")
        h_cov_reco = (TH2D*) h_cov_flux->Clone();
    else
        h_cov_reco = (TH2D*) h_cov_xsec->Clone();

    // Convert the Covariance Matrix-- switching from MC CV deviations to Data CV deviation
    _util.ConvertCovarianceUnits(h_cov_reco,
                            h_dataxsec_model.at(k_model_input),
                            h_dataxsec);

    // Now set the bin errors
    for (int bin = 1; bin < h_dataxsec->GetNbinsX()+1; bin++ ){    
        h_dataxsec->SetBinError(bin, std::sqrt(h_cov_reco->GetBinContent(bin, bin)));
    }

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.15);
    h_dataxsec->GetYaxis()->SetTitleOffset(1.7);
    h_dataxsec->Draw("E1,X0");
    
    // Draw the model data xsections
    for (unsigned int m = 0; m < models.size(); m++){

         if (m == k_model_input)
            continue;

        if (genmode == "flux"){
            h_dataxsec_model.at(k_model_FLUGG)->Draw("hist,same");
            break;
        }
        else {
            // Skip Tune 1
            if (m == k_model_tune1 || m  == k_model_FLUGG)
                continue;

            h_dataxsec_model.at(m)->Draw("hist,same");
        }

    }
    
    h_dataxsec->Draw("E1,same,X0");
    

    TLegend *leg = new TLegend(0.5, 0.5, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_dataxsec, "Data (Sys.)", "ep");
    
    if (genmode == "flux"){
        leg->AddEntry(h_dataxsec_model.at(k_model_FLUGG)    , "Data FLUGG", "l");
    }
    else {
        leg->AddEntry(h_dataxsec_model.at(k_model_mec)      , "Data 1.5 #times MEC", "l");
        leg->AddEntry(h_dataxsec_model.at(k_model_geniev3)  , "Data GENIE v3", "l");
        leg->AddEntry(h_dataxsec_model.at(k_model_nuwro)    , "Data NuWro", "l");
        // leg->AddEntry(h_dataxsec_model.at(k_model_tune1)    , "Data Tune 1", "l");
    }

    
    leg->Draw();

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.52, 0.92, 0.52, 0.92);

    if (genmode == "flux")
        c->Print(Form("plots/run%s/Models/%s/run%s_ModelDataComparison_%s_FLUGG.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));
    else
        c->Print(Form("plots/run%s/Models/%s/run%s_ModelDataComparison_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));
    
    delete c;

    fxsec->Close();

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareTotalDataCrossSections(){

    gStyle->SetOptStat(0);

    TH1D* h_temp;

    std::string genmode = "flux";
    // std::string genmode = "other";


    // Create a vector for the models
    std::vector<std::string> models = {
        "mec",
        "geniev3",
        "nuwro",
        "FLUGG",
        "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_mec,
        k_model_geniev3,
        k_model_nuwro,
        k_model_FLUGG,
        k_model_tune1,
        k_MODEL_MAX
    };

    std::vector<TH1D*> h_dataxsec_model(models.size());

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Data XSec
    h_temp  = (TH1D*)fxsec->Get(Form("CV/integrated/h_run%s_CV_0_integrated_data_xsec", _util.run_period));
    TH1D* h_dataxsec = (TH1D*)h_temp->Clone();

    // Set the error to be 21% total systematic uncertainty
    if (genmode == "flux")
        h_dataxsec->SetBinError(1,h_dataxsec->GetBinContent(1) * 0.216 );
    else
        h_dataxsec->SetBinError(1,h_dataxsec->GetBinContent(1) * std::sqrt(0.032*0.032 + 0.056*0.056) );
    
    // h_dataxsec->SetBinError(1, h_dataxsec->GetBinContent(1) * 0.21);
    h_dataxsec->SetLineColor(kBlack);
    h_dataxsec->SetMarkerStyle(20);
    h_dataxsec->SetMarkerSize(0.5);

    // X-Axis
    h_dataxsec->GetXaxis()->SetRangeUser(0.0,1.0); 
    h_dataxsec->GetXaxis()->SetLabelOffset(999);
    h_dataxsec->GetXaxis()->SetLabelSize(0);
    h_dataxsec->GetXaxis()->SetTickLength(0);
    
    // Y-Axis
    h_dataxsec->GetYaxis()->SetRangeUser(3.0, 8.0);
    h_dataxsec->GetYaxis()->CenterTitle();
    h_dataxsec->GetYaxis()->SetLabelSize(0.1);
    h_dataxsec->GetYaxis()->SetTitleSize(0.1);

    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
    
        // Data X Sec MEC
        h_temp  = (TH1D*)fxsec->Get(Form("%s/integrated/h_run%s_CV_0_integrated_data_xsec", models.at(m).c_str(), _util.run_period));
        h_dataxsec_model.at(m) = (TH1D*)h_temp->Clone();
        h_dataxsec_model.at(m)->SetLineWidth(2);

        // Set the line colours
        if (m == k_model_mec)      h_dataxsec_model.at(k_model_mec)      ->SetLineColor(kGreen+2);
        if (m == k_model_geniev3)  h_dataxsec_model.at(k_model_geniev3)  ->SetLineColor(kBlue+2);
        if (m == k_model_nuwro)    h_dataxsec_model.at(k_model_nuwro)    ->SetLineColor(kPink+1);
        if (m == k_model_FLUGG)    h_dataxsec_model.at(k_model_FLUGG)    ->SetLineColor(kViolet-1);
        if (m == k_model_tune1)    h_dataxsec_model.at(k_model_tune1)    ->SetLineColor(kOrange-1);
    
    }

    TCanvas *c = new TCanvas("c", "c", 150, 350);
    gPad->SetLeftMargin(0.3);
    h_dataxsec->Draw("E1,X0");
    
    // Draw the model data xsections
    for (unsigned int m = 0; m < models.size(); m++){
        
        if (genmode == "flux"){
            h_dataxsec_model.at(k_model_FLUGG)->Draw("hist,same");
            break;
        }
        else {
            // Skip Tune 1
            if (m == k_model_tune1 || m  == k_model_FLUGG)
                continue;

            h_dataxsec_model.at(m)->Draw("hist,same");
        }
    }
    
    h_dataxsec->Draw("E1,same,X0");
    h_dataxsec->SetMarkerSize(0.4);
    

    TLegend *leg = new TLegend(0.35, 0.70, 0.70, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_dataxsec, "Data (Sys.)", "ep");
    
    if (genmode == "flux"){
        leg->AddEntry(h_dataxsec_model.at(k_model_FLUGG)    , "Data FLUGG", "l");
    }
    else {
        leg->AddEntry(h_dataxsec_model.at(k_model_mec)      , "Data 1.5 #times MEC", "l");
        leg->AddEntry(h_dataxsec_model.at(k_model_geniev3)  , "Data GENIE v3", "l");
        leg->AddEntry(h_dataxsec_model.at(k_model_nuwro)    , "Data NuWro", "l");
        // leg->AddEntry(h_dataxsec_model.at(k_model_tune1)    , "Data Tune 1", "l");
    }
    leg->Draw();

    gStyle->SetLegendTextSize(0.06);

    if (genmode == "flux")
        c->Print(Form("plots/run%s/Models/Total/run%s_ModelDataComparison_FLUGG.pdf", _util.run_period, _util.run_period));
    else 
        c->Print(Form("plots/run%s/Models/Total/run%s_ModelDataComparison.pdf", _util.run_period, _util.run_period));
    
    delete c;

    fxsec->Close();

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareUnfoldedDataCrossSections(){

    gStyle->SetOptStat(0);

    // Create a vector for the models
    std::vector<std::string> models = {
        "Input",
        "mec",
        "geniev3",
        "nuwro",
        "FLUGG"
        // "tune1"
    };

    // enums for the models
    enum enum_models {
        k_model_input,
        k_model_mec,
        k_model_geniev3,
        k_model_nuwro,
        k_model_FLUGG,
        // k_model_tune1,
        k_MODEL_MAX
    };

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH2D* h_temp_2D;
    TH1D* h_temp;

    // Total Covariance Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_cov_tot_dataxsec_reco",_util.xsec_var));
    TH2D* h_cov_reco = (TH2D*)h_temp_2D->Clone();
    h_cov_reco->SetDirectory(0);

    // Flux Covariance Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_cov_flux_sys_dataxsec_reco",_util.xsec_var));
    TH2D* h_cov_reco_flux = (TH2D*)h_temp_2D->Clone();
    h_cov_reco_flux->SetDirectory(0);

    // Xsec Covariance Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_cov_xsec_sys_dataxsec_reco",_util.xsec_var));
    TH2D* h_cov_reco_xsec = (TH2D*)h_temp_2D->Clone();
    h_cov_reco_xsec->SetDirectory(0);

    fxsec->Close();

    // Now Get the Models
    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    std::vector<TH1D*> h_true_model(k_MODEL_MAX);
    std::vector<TH1D*> h_data_model(k_MODEL_MAX);
    std::vector<TH2D*> h_cov_diag(k_MODEL_MAX);
    std::vector<TH2D*> h_response_model(k_MODEL_MAX);
    std::vector<TH1D*> h_unf_model(k_MODEL_MAX);

    // Make the plot
    TCanvas *c = new TCanvas("c", "c", 500, 500);
    
    // Create the legend
    TLegend *leg = new TLegend(0.35, 0.50, 0.70, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){

        // Get true model xsec
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run%s_CV_0_%s_mc_xsec", models.at(m).c_str(),_util.vars.at(k_var_trueX).c_str(), _util.run_period, _util.vars.at(k_var_trueX).c_str()));
        h_true_model.at(m) = (TH1D*)h_temp->Clone();

        // Get data cross section extracted for model 
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run%s_CV_0_%s_data_xsec", models.at(m).c_str(), _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
        h_data_model.at(m) = (TH1D*)h_temp->Clone();

        // Get the response matrix for model
        h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run%s_%s_0_smearing", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), _util.run_period, models.at(m).c_str()));
        h_response_model.at(m) = (TH2D*)h_temp_2D->Clone();

        // Flip the response matrix 
        // Loop over rows
        for (int row=0; row<h_temp_2D->GetXaxis()->GetNbins()+2; row++) {

            for (int col=0; col<h_temp_2D->GetYaxis()->GetNbins()+2; col++){
                h_response_model.at(m)->SetBinContent(col, row, h_temp_2D->GetBinContent(row, col));          
            }
        }

        // Clone covariance matrix
        // Set diagonals of covariance matrix
        // if (m == k_model_FLUGG)
        //     h_cov_diag.at(m) = (TH2D*)h_cov_reco_xsec->Clone();
        // else if (m == k_model_input)
        //     h_cov_diag.at(m) = (TH2D*)h_cov_reco->Clone();
        // else
        //     h_cov_diag.at(m) = (TH2D*)h_cov_reco_xsec->Clone();
        
        h_cov_diag.at(m) = (TH2D*)h_cov_reco->Clone();
        
        // Initialise the WienerSVD class
        WienerSVD _wSVD;
        _wSVD.Initialise(_util);
        _wSVD.DoUnfolding(2, 0, h_true_model.at(m), h_data_model.at(m), h_response_model.at(m), h_cov_diag.at(m));

        h_unf_model.at(m) = (TH1D*)_wSVD.unf->Clone(Form("test_%s", models.at(m).c_str()));

        for (int bin = 1; bin < h_unf_model.at(m)->GetNbinsX()+1; bin++){
            double err = _wSVD.unfcov->GetBinContent(bin, bin);
            h_unf_model.at(m)->SetBinError(bin, std::sqrt(err));
        }
        
        // Set the line colours
        if (m == k_model_input){
            h_unf_model.at(m)->SetLineColor(kBlack);
            h_unf_model.at(m)->SetMarkerStyle(20);
            h_unf_model.at(m)->SetMarkerSize(0.5);
        }    
        
        if (m == k_model_mec)      h_unf_model.at(m)->SetLineColor(kGreen+2);
        if (m == k_model_geniev3)  h_unf_model.at(m)->SetLineColor(kBlue+2);
        if (m == k_model_nuwro)    h_unf_model.at(m)->SetLineColor(kPink+1);
        if (m == k_model_FLUGG)    h_unf_model.at(m)->SetLineColor(kViolet-1);
        // if (m == k_model_tune1)    h_unf_model.at(m)->SetLineColor(kOrange-1);
        
        
        h_unf_model.at(m)->Scale(1.0, "width");
        
        if (m == k_model_input){
            _util.IncreaseLabelSize( h_unf_model.at(m), c);
            gPad->SetLeftMargin(0.20);
            c->SetBottomMargin(0.15);
            h_unf_model.at(m)->SetTitle(_util.var_labels_xsec.at(k_var_trueX).c_str());
            h_unf_model.at(m)->Draw("E,same");
            h_unf_model.at(m)->GetYaxis()->SetTitleOffset(1.4);
        }
        else {
            h_unf_model.at(m)->Draw("hist,same");
        }

        h_unf_model.at(m)->SetLineWidth(2);

        h_unf_model.at(m)->SetMinimum(0);

        
        if (m == k_model_input) {
            models.at(m) = "CV";
            leg->AddEntry(h_unf_model.at(m), "Data (Stat. + Sys.)", "ep");
        }
        if (m == k_model_mec)      leg->AddEntry(h_unf_model.at(m), "Data 1.5 #times MEC", "l");
        if (m == k_model_geniev3)  leg->AddEntry(h_unf_model.at(m), "Data GENIE v3", "l");
        if (m == k_model_nuwro)    leg->AddEntry(h_unf_model.at(m), "Data NuWro", "l");
        if (m == k_model_FLUGG)    leg->AddEntry(h_unf_model.at(m), "Data FLUGG", "l");
        // if (m == k_model_tune1)    leg->AddEntry(h_unf_model.at(m), "Data Tune 1", "l");
        

        delete _wSVD.smear;
        delete _wSVD.wiener;
        delete _wSVD.unfcov;
        delete _wSVD.unf;
        delete _wSVD.diff;
        delete _wSVD.bias;
        delete _wSVD.bias2;
        delete _wSVD.fracError;
        delete _wSVD.absError;
        delete _wSVD.MSE;
        delete _wSVD.MSE2;
    }

    h_unf_model.at(k_model_input)->Draw("E,same");
    gStyle->SetLegendTextSize(0.04);
    leg->Draw();

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.52, 0.92, 0.52, 0.92);

    c->Print(Form("plots/run%s/Models/%s/run%s_UnfoldedDataComparison_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));
    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::SaveResponseMatrix(){

    TFile *fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    _util.CreateDirectory("Response/");

    TH2D* h_temp_2D;

    std::vector<std::string> variables = {"elec_E", "elec_ang", "elec_cang"};

    // For smearing
    // std::vector<std::string> var_names = {"reco_el_E", "reco_el_ang", "reco_el_cang"};

    // For response
    std::vector<std::string> var_names = {"true_el_E", "true_el_ang", "true_el_cang"};

    // Vector for storing the response matrix
    std::vector<TH2D*> h_response(variables.size());
    std::vector<TH2D*> h_response_index(variables.size());

    for (unsigned int m = 0; m < variables.size(); m++){
        
        // Get the response matrix
        h_temp_2D  = (TH2D*)fxsec->Get(Form("CV/%s/h_run1_CV_0_smearing", var_names.at(m).c_str()));
        h_response.at(m) = (TH2D*)h_temp_2D->Clone();
    
        // Now we got the histogram lets save it!
        TCanvas * c = new TCanvas("c", "c", 500, 500);
        c->SetTopMargin(0.11);

        h_response.at(m)->SetStats(kFALSE);

        _util.IncreaseLabelSize(h_response.at(m), c);

        gStyle->SetPalette(kBlueGreenYellow);
        gStyle->SetPaintTextFormat("4.2f");
        h_response.at(m)->SetMarkerSize(0.4);
        h_response.at(m)->SetMarkerColor(kRed+2);

        h_response.at(m)->Draw("colz");

        // Draw the run period on the plot
        // _util.Draw_Run_Period(c, 0.76, 0.915, 0.76, 0.915);

        c->Print(Form("plots/run%s/Response/Response_run%s_%s.pdf", _util.run_period, _util.run_period, variables.at(m).c_str()));
        delete c;
    
    }

    // Now convert the histogram to bin indexes to read the response matrix more clearly
    for (unsigned int m = 0; m < variables.size(); m++){
        
        // Get the response matrix
        h_response_index.at(m) = new TH2D("", ";True Bin i; Reconstructed Bin j", h_response.at(m)->GetNbinsX(), 1, h_response.at(m)->GetNbinsX()+1, h_response.at(m)->GetNbinsY(), 1, h_response.at(m)->GetNbinsY()+1);
    
        // Set the bin values
        for (int x = 1; x < h_response.at(m)->GetNbinsY()+1; x++){
            for (int y = 1; y < h_response.at(m)->GetNbinsX()+1; y++){
                 h_response_index.at(m)->SetBinContent(x,y,h_response.at(m)->GetBinContent(x,y));
            }
        }

        if (variables.at(m) == "elec_E"){
            h_response_index.at(m)->SetTitle("E_{e}");
        }
        if (variables.at(m) == "elec_ang"){
            h_response_index.at(m)->SetTitle("#beta_{e#lower[-0.5]{-} + e^{+}}");
        }
        if (variables.at(m) == "elec_cang"){
            h_response_index.at(m)->SetTitle("cos#beta_{e}");
        }

        // Write the histograms to file
        TFile *f_result = TFile::Open("files/xsec_result_run1_paper.root", "UPDATE");
        
        if (variables.at(m) == "elec_E"){
            h_response_index.at(m)->Write("response_energy",TObject::kOverwrite);
        }
        if (variables.at(m) == "elec_cang"){
            h_response_index.at(m)->Write("response_angle",TObject::kOverwrite);
        }
    
        f_result->Close();



        // Now we got the histogram lets save it!
        TCanvas * c = new TCanvas("c", "c", 500, 500);
        c->SetTopMargin(0.11);

        h_response_index.at(m)->SetStats(kFALSE);

        _util.IncreaseLabelSize(h_response_index.at(m), c);

        gStyle->SetPalette(kBlueGreenYellow);
        gStyle->SetPaintTextFormat("4.2f");
        h_response_index.at(m)->SetMarkerSize(1.0);
        h_response_index.at(m)->SetMarkerColor(kRed+2);

        h_response_index.at(m)->Draw("colz,text00");
        h_response_index.at(m)->GetXaxis()->CenterLabels(kTRUE);
        h_response_index.at(m)->GetYaxis()->CenterLabels(kTRUE);
        h_response_index.at(m)->GetXaxis()->SetNdivisions(h_response.at(m)->GetNbinsX(), 0, 0, kFALSE);
        h_response_index.at(m)->GetYaxis()->SetNdivisions(h_response.at(m)->GetNbinsY(), 0, 0, kFALSE);


        // Draw the run period on the plot 
        // _util.Draw_Run_Period(c, 0.76, 0.915, 0.76, 0.915);
        _util.Draw_ubooneSim(c, 0.33, 0.925, 0.33, 0.905);

        c->Print(Form("plots/run%s/Response/Response_run%s_%s_index.pdf", _util.run_period, _util.run_period, variables.at(m).c_str()));
        delete c;
    
    }

    // MC XSec
    TH1D* h_temp;
    h_temp  = (TH1D*)fxsec->Get(Form("geniev3/%s/h_run1_CV_0_%s_mc_xsec", _util.vars.at(k_var_trueX).c_str(), _util.vars.at(k_var_trueX).c_str()));
    TH1D* truexsec = (TH1D*)h_temp->Clone();
    truexsec->Scale(1.0, "width");
    truexsec->SetDirectory(0);
    fxsec->Close();

    TFile *f_result = TFile::Open("files/xsec_result_run1_paper.root", "UPDATE");
    if (std::string(_util.xsec_var) == "elec_E"){
        truexsec->Write("mc_xsec_true_genie_v3_0_6_energy",TObject::kOverwrite);
    }
    else {
        truexsec->Write("mc_xsec_true_genie_v3_0_6_angle",TObject::kOverwrite);
    }
    f_result->Close();


   

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CheckPi0Coverage(){

    gStyle->SetOptStat(0);

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH2D* h_temp_2D;
    TH1D* h_temp;

    // The covariance matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/er/h_cov_genie_multi_mcxsec_reco",_util.xsec_var));
    TH2D* h_cov = (TH2D*)h_temp_2D->Clone();
    h_cov->SetDirectory(0);

    fxsec->Close();

    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // MC R CV
    h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", "CV", _util.vars.at(k_var_recoX).c_str(), _util.vars.at(k_var_recoX).c_str()));
    TH1D* h_CV = (TH1D*)h_temp->Clone();
    h_CV->Scale(1.0, "width");

    // MC R CV
    h_temp  = (TH1D*)fxsec->Get(Form("nopi0tune/%s/h_run1_CV_0_%s_mc_xsec", _util.vars.at(k_var_recoX).c_str(), _util.vars.at(k_var_recoX).c_str()));
    TH1D* h_nopi0 = (TH1D*)h_temp->Clone();
    h_nopi0->Scale(1.0, "width");


    // Set the bin error of the cv to be for the genie multisim systematic uncertainty
    for (int bin = 1; bin < h_CV->GetNbinsX(); bin++){
        h_CV->SetBinError(bin, std::sqrt(h_cov->GetBinContent(bin, bin)));
    }

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.15);
    h_CV->GetYaxis()->SetTitleOffset(1.7);

    // Set the line properties
    h_CV->SetLineColor(kBlack);
    h_nopi0->SetLineColor(kPink+1);

    h_CV->SetLineWidth(2);
    h_nopi0->SetLineWidth(2);

    // Draw the CV
    h_CV->Draw("hist,E");
    h_nopi0->Draw("hist,same");
    h_CV->Draw("hist,E,same");
    
    TLegend *leg = new TLegend(0.4, 0.5, 0.6, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(h_CV, "MC CV (Genie All Sys.)", "el");
    leg->AddEntry(h_nopi0, "MC no #pi^{0} Tune", "l");
    leg->Draw();

    c->Print(Form("plots/run%s/Systematics/pi0/%s/run%s_pi0tune_sys_coverage_%s.pdf", _util.run_period, _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
    delete c;


}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareMCC8Result(){

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c", 200, 350);
    gPad->SetLeftMargin(0.3);
    
    // Data X-Sec with Stat Only
    TH1D* h_data = new TH1D("h_data", ";;#nu_{e} + #bar{#nu}_{e} CC Cross Section [cm^{2} / nucleon]", 2, 0, 2);
    TH1D* h_data_mcc9 = new TH1D("h_data2", ";;#nu_{e} + #bar{#nu}_{e} CC Cross Section [cm^{2} / nucleon]", 2, 0, 2);
    
    // X-Axis
    h_data->GetXaxis()->SetRangeUser(0.0,2.0); 
    // h_data->GetXaxis()->SetLabelOffset(999);
    h_data->GetXaxis()->SetLabelSize(0.08);
    h_data->GetXaxis()->SetTickLength(0);
    
    // Y-Axis
    h_data->GetYaxis()->SetRangeUser(0.22e-38,1.22e-38);
    h_data->GetYaxis()->CenterTitle();
    h_data->GetYaxis()->SetLabelSize(0.1);
    h_data->GetYaxis()->SetTitleSize(0.08);
   
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.7);
    h_data->SetLineColor(kBlack);

    h_data_mcc9->SetMarkerStyle(23);
    h_data_mcc9->SetMarkerSize(0.7);
    h_data_mcc9->SetLineColor(kBlack);


    // Fill it
    h_data->Fill("MCC8", 6.8426915e-39); // new in FV flux
    h_data->Fill("MCC9", 0.0); // new in FV flux
    h_data->SetBinError(1, 6.8426915e-39 * 0.40); // new in FV flux
    h_data->SetBinError(2, 0.0); // new in FV flux
    h_data->Draw("E1,X0");
    
    // Stat Error
    h_data_mcc9->Fill("MCC8", 0.0); // new in FV flux
    h_data_mcc9->Fill("MCC9", 6.62447e-39); // new in FV flux
    h_data_mcc9->SetBinError(1, 0.0); // new in FV flux
    h_data_mcc9->SetBinError(2, 6.62447e-39 * 0.244); // new in FV flux
    h_data_mcc9->Draw("E1,X0,same");

    // Statistical band
    TH1D * h_data_stat = (TH1D*) h_data->Clone();
    // h_data_stat->SetBinError(1, 0.144e-38);
    h_data_stat->SetBinError(1, 6.8426915e-39 * 0.22 ); // new in FV flux
    h_data_stat->SetLineColor(kBlack);
    h_data_stat->Draw("E1,X0,same");

    // Statistical band
    TH1D * h_data_stat_mcc9 = (TH1D*) h_data_mcc9->Clone();
    // h_data_stat->SetBinError(1, 0.144e-38);
    h_data_stat_mcc9->SetBinError(2, 6.62447e-39 * 0.0987 ); // new in FV flux
    h_data_stat_mcc9->SetLineColor(kBlack);
    h_data_stat_mcc9->Draw("E1,X0,same");


    // Genie v12.2.2 nue + nuebar
    TH1D* h_genie_v2_nue_nuebar = new TH1D("h_genie_v2", "", 1, 0.0, 2.0);
    // h_genie_v2_nue_nuebar->Fill(0.5,7.19925e-39 );
    h_genie_v2_nue_nuebar->Fill(0.5,7.3125100e-39 ); // with FV flux
    h_genie_v2_nue_nuebar->SetLineColor(kViolet-5);
    h_genie_v2_nue_nuebar->SetLineWidth(3); 
    h_genie_v2_nue_nuebar->SetLineStyle(7);
    h_genie_v2_nue_nuebar->Draw("hist,same");

    // Genie v3 nue + nuebar
    TH1D* h_genie_v3_nue_nuebar = new TH1D("h_genie_v3", "", 1, 0.0, 2.0);
    // h_genie_v3_nue_nuebar->Fill(0.5,5.5228738e-39 );
    h_genie_v3_nue_nuebar->Fill(0.5,5.5711475e-39 ); // with FV flux
    h_genie_v3_nue_nuebar->SetLineColor(kBlue+2);
    h_genie_v3_nue_nuebar->SetLineWidth(3);
    h_genie_v3_nue_nuebar->SetLineStyle(8);
    h_genie_v3_nue_nuebar->Draw("hist,same");

    // NuWro nue + nuebar
    TH1D* h_genie_NuWro_nue_nuebar = new TH1D("h_nuwro_v2", "", 1, 0.0, 2.0);
    // h_genie_NuWro_nue_nuebar->Fill(0.5,3.8158e-39 );
    h_genie_NuWro_nue_nuebar->Fill(0.5,5.9205940e-39 ); // with FV flux
    h_genie_NuWro_nue_nuebar->SetLineColor(kRed+2);
    h_genie_NuWro_nue_nuebar->SetLineWidth(3);
    h_genie_NuWro_nue_nuebar->SetLineStyle(1);
    h_genie_NuWro_nue_nuebar->Draw("hist,same");


    h_data->Draw("E1,X0,same");
    h_data_stat->Draw("E1,X0,same");


    // Draw the Legend
    TLegend *leg = new TLegend(0.35, 0.70, 0.70, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data, "Data (stat. + sys.)",      "ep");
    leg->AddEntry(h_data_mcc9, "Data (stat. + sys.)",      "ep");
    leg->AddEntry(h_genie_v2_nue_nuebar,   "GENIE v2.12.2",    "l");
    leg->AddEntry(h_genie_v3_nue_nuebar,   "GENIE v3.0.6",    "l");
    leg->AddEntry(h_genie_NuWro_nue_nuebar,   "NuWro v19.02.1",    "l");
    
    leg->Draw();

    gStyle->SetLegendTextSize(0.06);
   
    TLatex *t = new TLatex(.34, .145, "#splitline{MicroBooNE NuMI}{Data}");
    t->SetTextColor(kBlack);
    t->SetNDC();
    t->SetTextSize(2.0/30.);
    t->SetTextAlign(11);
    // t->Draw();


    c->Print("plots/mcc8_mcc9_nuexsec_generator_plot.pdf");

}
// -----------------------------------------------------------------------------
void UtilityPlotter::ForwardFoldedGeneratorComparison(){

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
        "geniev3",
        "geniev2gen",
        "nuwrogen"
    };

    // enums for the models
    enum enum_models {
        k_model_CV,
        k_model_geniev3,
        k_model_geniev2gen,
        k_model_nuwrogen,
        k_MODEL_MAX
    };

    // Create the vector of histograms
    std::vector<TH2D*> h_response_model(models.size());
    std::vector<TH1D*> h_mcxsec_true_model(models.size());
    std::vector<TH1D*> h_mcxsec_reco_model(models.size());


    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
        
        // Response Matrix -- only get the CV here
        h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run1_%s_0_smearing", models.at(k_model_CV).c_str(), _util.vars.at(k_var_trueX).c_str(), models.at(k_model_CV).c_str()));
        if (h_temp_2D == NULL) std::cout <<"Help!" << m << std::endl;
        h_response_model.at(m) = (TH2D*)h_temp_2D->Clone();

        // MC xsec in True
        if (m == k_model_geniev2gen || m == k_model_nuwrogen){
            h_temp  = (TH1D*)fxsec->Get(Form("%s/h_%s", models.at(m).c_str(), _util.xsec_var));
            std::cout << Form("%s/h_%s", models.at(m).c_str(), _util.xsec_var) << std::endl;
        }
        else
            h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), _util.vars.at(k_var_trueX).c_str()));
        
        h_mcxsec_true_model.at(m) = (TH1D*)h_temp->Clone();
        h_mcxsec_reco_model.at(m) = (TH1D*)h_dataxsec->Clone();
       
        _util.MatrixMultiply(h_mcxsec_true_model.at(m), h_mcxsec_reco_model.at(m), h_response_model.at(k_model_CV), "true_reco", true);

        h_mcxsec_reco_model.at(m)->SetLineWidth(3);
    }

    // Set the line colours
    h_mcxsec_reco_model.at(k_model_CV)        ->SetLineColor(kRed+2);
    h_mcxsec_reco_model.at(k_model_geniev3)   ->SetLineColor(kGreen+2);
    h_mcxsec_reco_model.at(k_model_geniev2gen)->SetLineColor(kOrange-1);
    h_mcxsec_reco_model.at(k_model_nuwrogen)  ->SetLineColor(kPink+1);

    // Set the line styles
    h_mcxsec_reco_model.at(k_model_CV)         ->SetLineStyle(1);
    h_mcxsec_reco_model.at(k_model_geniev3)    ->SetLineStyle(2);
    h_mcxsec_reco_model.at(k_model_geniev2gen) ->SetLineStyle(3);
    h_mcxsec_reco_model.at(k_model_nuwrogen)   ->SetLineStyle(4);
    

    // Now lets plot
    TCanvas *c = new TCanvas("c", "c", 500, 500);
    c->SetLeftMargin(0.2);
    c->SetBottomMargin(0.15);
    h_dataxsec->GetYaxis()->SetTitleOffset(1.7);
    h_dataxsec->SetMinimum(0);
    h_dataxsec->Draw("E1,X0,same");

    if ( std::string(_util.xsec_var) == "elec_cang"){
        h_dataxsec->GetYaxis()->SetRangeUser(0, 10);
        h_dataxsec->GetXaxis()->SetLabelSize(0.03);
    }
    
    h_mcxsec_reco_model.at(k_model_CV)->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_geniev3)->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_geniev2gen)->Draw("hist,same");
    h_mcxsec_reco_model.at(k_model_nuwrogen)->Draw("hist,same");
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
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_CV),  Form("MC #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");
    
    std::cout << "Genie v3" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_geniev3), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_geniev3),  Form("GENIE v3.0.6 #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");

    std::cout << "Geniev2" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_geniev2gen), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_geniev2gen),  Form("GENIE v2.12.2 #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");

    std::cout << "NuWro" << std::endl;
    _util.CalcChiSquared(h_mcxsec_reco_model.at(k_model_nuwrogen), h_dataxsec, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_reco_model.at(k_model_nuwrogen),  Form("NuWro v19.02.02 #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");


    gStyle->SetLegendTextSize(0.03);
    
    leg->Draw();

    // Draw the run period on the plot
    // _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.52, 0.92, 0.52, 0.92);
    
    c->Print(Form("plots/run%s/Models/%s/run%s_DataGeneratorComparison_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));

    fxsec->Close();

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareGeneratorTotalCrossSec(){

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c", 150, 350);
    gPad->SetLeftMargin(0.3);

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH1D* h_temp;

    // Data Xsec stat
    h_temp  = (TH1D*)fxsec->Get("total/h_data_xsec_stat_reco");
    TH1D* h_data_stat = (TH1D*)h_temp->Clone();
    h_data_stat->SetDirectory(0);

    // Data Xsec stat + sys
    h_temp  = (TH1D*)fxsec->Get("total/h_data_xsec_stat_sys_reco");
    TH1D* h_data = (TH1D*)h_temp->Clone();
    h_data->SetDirectory(0);

    fxsec->Close();

    // Now Get the Models
    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    // Now get some other models
    // Create a vector for the models
    // Create a vector for the models
    std::vector<std::string> models = {
        "CV",
        "geniev3",
        "geniev2gen",
        "nuwrogen"
    };

    // enums for the models
    enum enum_models {
        k_model_CV,
        k_model_geniev3,
        k_model_geniev2gen,
        k_model_nuwrogen,
        k_MODEL_MAX
    };

    std::vector<TH1D*> h_model_xsec(k_MODEL_MAX);
    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){

        // Get true tune1 xsec
        if (m == k_model_geniev2gen || m == k_model_nuwrogen)
            h_temp  = (TH1D*)fxsec->Get(Form("%s/h_elec_tot", models.at(m).c_str()));
        else
            h_temp  = (TH1D*)fxsec->Get(Form("%s/integrated/h_run%s_CV_0_integrated_mc_xsec", models.at(m).c_str(), _util.run_period));
        
        h_model_xsec.at(m) = (TH1D*)h_temp->Clone();
        h_model_xsec.at(m)->SetLineWidth(3);
       
    }


    h_model_xsec.at(k_model_CV)->SetLineColor(kRed+2);
    h_model_xsec.at(k_model_geniev3)->SetLineColor(kGreen+2);
    h_model_xsec.at(k_model_geniev2gen)->SetLineColor(kOrange-1);
    h_model_xsec.at(k_model_nuwrogen)->SetLineColor(kPink+1);

    // Set the line styles
    h_model_xsec.at(k_model_CV)       ->SetLineStyle(1);
    h_model_xsec.at(k_model_geniev3)  ->SetLineStyle(2);
    h_model_xsec.at(k_model_geniev2gen)    ->SetLineStyle(3);
    h_model_xsec.at(k_model_nuwrogen)    ->SetLineStyle(4);


    // X-Axis
    h_data->GetXaxis()->SetRangeUser(0.0,1.0); 
    h_data->GetXaxis()->SetLabelOffset(999);
    h_data->GetXaxis()->SetLabelSize(0);
    h_data->GetXaxis()->SetTickLength(0);
    
    // Y-Axis
    h_data->GetYaxis()->SetRangeUser(3.0, 8.0);
    h_data->GetYaxis()->CenterTitle();
    h_data->GetYaxis()->SetLabelSize(0.1);
    h_data->GetYaxis()->SetTitleSize(0.1);
   
    // h_data->SetLineWidth(2);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.7);
    h_data->SetLineColor(kBlack);
    h_data->Draw("E1,X0");

    // Statistical band
    h_data_stat->SetLineColor(kBlack);
    h_data_stat->Draw("E1,X0,same");

    // Draw the models
    for (unsigned int m = 0; m < models.size(); m++){
        h_model_xsec.at(m)->Draw("hist,same");
    }


    h_data->Draw("E1,X0,same");
    h_data_stat->Draw("E1,X0,same");


    // Draw the Legend
    TLegend *leg = new TLegend(0.35, 0.70, 0.70, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    h_data->SetMarkerSize(0.4);
    leg->AddEntry(h_data, "Data (stat. + sys.)",        "ep");
    leg->AddEntry( h_model_xsec.at(k_model_CV),         "MC", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_geniev3),    "GENIE v3.0.6", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_geniev2gen), "GENIE v2.12.2", "lf");
    leg->AddEntry( h_model_xsec.at(k_model_nuwrogen),      "NuWro v19.02.02", "lf");
    
    leg->Draw();

    gStyle->SetLegendTextSize(0.06);

    double Data_POT; 

    // Set the scale factors
    if (strcmp(_util.run_period, "1") == 0){
        Data_POT = _util.config_v.at(_util.k_Run1_Data_POT); // Define this variable here for easier reading
    }
    else if (strcmp(_util.run_period, "3") == 0){
        Data_POT = _util.config_v.at(_util.k_Run3_Data_POT); // Define this variable here for easier reading
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
        exit(1);
    }

    Data_POT = Data_POT / 1.0e20;

    TLatex *t = new TLatex(.34, .145, Form("#splitline{MicroBooNE NuMI}{Data %2.1f#times10^{20} POT}", Data_POT));
    t->SetTextColor(kBlack);
    t->SetNDC();
    t->SetTextSize(2.0/30.);
    t->SetTextAlign(11);
    t->Draw();


    c->Print(Form("plots/run%s/Models/Total/run%s_GeneratorTotalCrossSectionComparison.pdf", _util.run_period, _util.run_period));

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareGeneratorUnfoldedModels(){

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
        "geniev3",
        "geniev2gen",
        "nuwro"
    };

    // enums for the models
    enum enum_models {
        k_model_CV,
        k_model_geniev3,
        k_model_geniev2gen,
        k_model_nuwro,
        k_MODEL_MAX
    };

    std::vector<TH1D*> h_mcxsec_true_model(models.size());
    std::vector<TH1D*> h_mcxsec_true_model_smear(models.size());

    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){
        // MC Xsec True
        
        if (m == k_model_geniev2gen)
            h_temp  = (TH1D*)fxsec->Get(Form("%s/h_%s", models.at(m).c_str(), _util.xsec_var));
        else
            h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run1_CV_0_%s_mc_xsec", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), _util.vars.at(k_var_trueX).c_str()));
        h_mcxsec_true_model.at(m)         = (TH1D*)h_temp->Clone();
        h_mcxsec_true_model_smear.at(m)   = (TH1D*)h_temp->Clone();

        _util.MatrixMultiply(h_mcxsec_true_model.at(m), h_mcxsec_true_model_smear.at(m), h_ac, "reco_true",false);

    }
    
   
    TLegend *leg = new TLegend(0.35, 0.6, 0.7, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(unf, "Data (Stat. + Sys.)", "ep");
    gStyle->SetLegendTextSize(0.03);
    
    // Now calculate the chi-squared
    double chi, pval;
    int ndof;

    std::cout << "CV" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_CV), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_CV),   Form("GENIE v3.0.6 (#muB tune) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");
    
    std::cout << "Genie v3" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_geniev3), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_geniev3),   Form("GENIE v3.0.6 #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");

    std::cout << "Genie v2" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_geniev2gen), unf, h_cov, chi, ndof, pval);
    // leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_geniev2gen),   Form("GENIE v2.12.2 #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");

    std::cout << "NuWro" << std::endl;
    _util.CalcChiSquared(h_mcxsec_true_model_smear.at(k_model_nuwro), unf, h_cov, chi, ndof, pval);
    leg->AddEntry(h_mcxsec_true_model_smear.at(k_model_nuwro),   Form("NuWro v19.02.2 #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "lf");

    // Scale the histograms by bin width 
    for (unsigned int m = 0; m < models.size(); m++){
        h_mcxsec_true_model_smear.at(m)->Scale(1.0, "width");
        h_mcxsec_true_model_smear.at(m)->SetLineWidth(5);
    }
    unf->Scale(1.0, "width");

    // Make the plot
    TCanvas *c = new TCanvas("c", "c", 500, 500);
    _util.IncreaseLabelSize( h_mcxsec_true_model_smear.at(k_model_CV), c);

    if ( std::string(_util.xsec_var) == "elec_cang"){
        h_mcxsec_true_model_smear.at(k_model_CV)->GetXaxis()->SetLabelSize(0.03);
    }

    gPad->SetLeftMargin(0.20);
    c->SetBottomMargin(0.15);
    h_mcxsec_true_model_smear.at(k_model_CV)->SetTitle(_util.var_labels_xsec.at(2).c_str());
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

    // Set the line styles
    h_mcxsec_true_model_smear.at(k_model_CV)         ->SetLineStyle(1);
    h_mcxsec_true_model_smear.at(k_model_geniev3)    ->SetLineStyle(2);
    h_mcxsec_true_model_smear.at(k_model_geniev2gen) ->SetLineStyle(3);
    h_mcxsec_true_model_smear.at(k_model_nuwro)      ->SetLineStyle(4);

    h_mcxsec_true_model_smear.at(k_model_CV)->SetMinimum(0.0);
    h_mcxsec_true_model_smear.at(k_model_CV)->Draw("hist");

    h_mcxsec_true_model_smear.at(k_model_geniev3)->SetLineColor(kGreen+2);
    h_mcxsec_true_model_smear.at(k_model_geniev3)->Draw("hist,same" );

    h_mcxsec_true_model_smear.at(k_model_geniev2gen)->SetLineColor(kOrange-1);
    // h_mcxsec_true_model_smear.at(k_model_geniev2gen)->Draw("hist,same" );

    h_mcxsec_true_model_smear.at(k_model_nuwro)->SetLineColor(kPink+1);
    h_mcxsec_true_model_smear.at(k_model_nuwro)->Draw("hist,same" );
    
    unf->Draw("E1,X0,same");

    // Draw the run period on the plot
    // _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.52, 0.92, 0.52, 0.92);
    

    leg->Draw();
    
    c->Print(Form("plots/run%s/Models/%s/run%s_DataModelGeneratorUnfoldedComparison_%s.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));

    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareGeneratorPi0(){

    // Load in the root file
    TFile *f_mc, *f_mc_nuwro, *f_mc_tune1;
    TTree *mc_tree, *mc_tree_nuwro, *mc_tree_tune1;

    // Get the TTree
    _util.GetFile(f_mc, "../ntuples/neutrinoselection_filt_run1_overlay_newtune.root"); // Get the run 1 MC file
    _util.GetTree(f_mc, mc_tree, "nuselection/NeutrinoSelectionFilter");

    SliceContainer SC;
    SC.Initialise(mc_tree, _util.k_mc, _util);

    SelectionCuts _scuts;
    _scuts.Initalise(_util);

    
    // 1D pi0 momentum
    TH1D *h_pi0_e_genie = new TH1D("h_true_pi0_e_genie", "; #pi^{0} Energy [GeV]; Entries", 15, 0, 1.5);
    TH1D *h_pi0_e_nuwro = new TH1D("h_true_pi0_e_nuwro", "; #pi^{0} Energy [GeV]; Entries", 15, 0, 1.5);
    TH1D *h_pi0_e_tune1 = new TH1D("h_true_pi0_e_tune1", "; #pi^{0} Energy [GeV]; Entries", 15, 0, 1.5);

    TH1D *h_pi0_e_genie_pass = new TH1D("h_true_pi0_e_genie_pass", "; #pi^{0} Energy [GeV]; Entries", 5, 0.0, 1.5);
    TH1D *h_pi0_e_nuwro_pass = new TH1D("h_true_pi0_e_nuwro_pass", "; #pi^{0} Energy [GeV]; Entries", 5, 0.0, 1.5);
    TH1D *h_pi0_e_tune1_pass = new TH1D("h_true_pi0_e_tune1_pass", "; #pi^{0} Energy [GeV]; Entries", 5, 0.0, 1.5);
    
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

        // Classify the event -- sets variable in the slice contianer
        SC.SliceClassifier(_util.k_mc);      // Classification of the event

        // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
        SC.SetThresholdEvent(_util.k_mc);

        // If the backtracked pdg of the leading shower is not an electron then alter classification
        SC.SetNonLdgShrEvent(_util.k_mc);
        
        SC.SliceInteractionType(_util.k_mc); // Genie interaction type
        SC.ParticleClassifier(_util.k_mc);   // The truth matched particle type of the leading shower
        SC.Pi0Classifier(_util.k_mc); 

        // Set derived variables in the slice container
        SC.SetSignal();                // Set the event as either signal or other
        SC.SetFakeData();              // Set the classifcation as data if fake data mode
        SC.SetTrueElectronThetaPhi();  // Set the true electron theta and phi variables
        SC.SetNuMIAngularVariables();  // Set the NuMI angular variables
        SC.CalibrateShowerEnergy();    // Divide the shower energy by 0.83 so it is done in one place

        bool is_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 

        double weight = _util.GetCVWeight(_util.k_mc, 1.0, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv, SC.interaction, SC.elec_e); // Turn off the genie tune
        
        // Pi0 Energy
        if (is_in_fv && SC.npi0>0 && SC.ccnc == 1) 
            h_pi0_e_genie->Fill(SC.pi0_e, weight);


        // Require Passed cuts
        if (_scuts.swtrig(SC, _util.k_mc) &&
        _scuts.slice_id(SC)            &&
        _scuts.e_candidate(SC)         &&
        _scuts.in_fv(SC)               &&
        _scuts.contained_frac(SC)      &&
        _scuts.topo_score(SC)          &&
        _scuts.shr_cosmic_IP(SC)       &&
        _scuts.shower_score(SC)        &&
        _scuts.shr_hitratio(SC)        &&
        _scuts.shr_moliere_avg(SC)     &&
        _scuts.shr_dist_dEdx_max(SC)   &&
        _scuts.dEdx_max_no_tracks(SC)){
            if (is_in_fv && SC.npi0>0 && SC.ccnc == 1) 
                h_pi0_e_genie_pass->Fill(SC.pi0_e, weight);
        }


    }

    // Get the TTree
    _util.GetFile(f_mc_nuwro, "../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_nuwro.root"); // Get the run 1 MC file
    _util.GetTree(f_mc_nuwro, mc_tree_nuwro, "nuselection/NeutrinoSelectionFilter");

    SliceContainer SC_nuwro;
    SC_nuwro.Initialise(mc_tree_nuwro, _util.k_mc, _util);

    int mc_tree_total_entries_nuwro = mc_tree_nuwro->GetEntries();
    std::cout << "Total MC Events:         " << mc_tree_total_entries_nuwro << std::endl;

    // Event loop
    for (int ievent = 0; ievent < mc_tree_total_entries_nuwro; ievent++){

        // See if we want to process all the events
        if (_util.num_events > 0){
            if (ievent >= _util.num_events) break;
        }

        // Alert the user
        if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
    
        // Get the entry in the tree
        mc_tree_nuwro->GetEntry(ievent); 

        // Classify the event -- sets variable in the slice contianer
        SC_nuwro.SliceClassifier(_util.k_mc);      // Classification of the event

        // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
        SC_nuwro.SetThresholdEvent(_util.k_mc);

        // If the backtracked pdg of the leading shower is not an electron then alter classification
        SC_nuwro.SetNonLdgShrEvent(_util.k_mc);
        
        SC_nuwro.SliceInteractionType(_util.k_mc); // Genie interaction type
        SC_nuwro.ParticleClassifier(_util.k_mc);   // The truth matched particle type of the leading shower
        SC_nuwro.Pi0Classifier(_util.k_mc); 

        // Set derived variables in the slice container
        SC_nuwro.SetSignal();                // Set the event as either signal or other
        SC_nuwro.SetFakeData();              // Set the classifcation as data if fake data mode
        SC_nuwro.SetTrueElectronThetaPhi();  // Set the true electron theta and phi variables
        SC_nuwro.SetNuMIAngularVariables();  // Set the NuMI angular variables
        SC_nuwro.CalibrateShowerEnergy();    // Divide the shower energy by 0.83 so it is done in one place

        SC_nuwro.SetPPFXCVWeight();

        bool is_in_fv = _util.in_fv(SC_nuwro.true_nu_vtx_sce_x, SC_nuwro.true_nu_vtx_sce_y, SC_nuwro.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 

        double weight = _util.GetCVWeight(_util.k_mc, SC_nuwro.weightSplineTimesTune, SC_nuwro.ppfx_cv, SC_nuwro.nu_e, SC_nuwro.nu_pdg, is_in_fv, SC_nuwro.interaction, SC_nuwro.elec_e);
        
        // Pi0 Energy
        if (is_in_fv && SC_nuwro.npi0>0 && SC_nuwro.ccnc == 1) 
            h_pi0_e_nuwro->Fill(SC_nuwro.pi0_e, weight);

        // Require Passed cuts
        if (_scuts.swtrig(SC_nuwro, _util.k_mc) &&
        _scuts.slice_id(SC_nuwro)            &&
        _scuts.e_candidate(SC_nuwro)         &&
        _scuts.in_fv(SC_nuwro)               &&
        _scuts.contained_frac(SC_nuwro)      &&
        _scuts.topo_score(SC_nuwro)          &&
        _scuts.shr_cosmic_IP(SC_nuwro)       &&
        _scuts.shower_score(SC_nuwro)        &&
        _scuts.shr_hitratio(SC_nuwro)        &&
        _scuts.shr_moliere_avg(SC_nuwro)     &&
        _scuts.shr_dist_dEdx_max(SC_nuwro)   &&
        _scuts.dEdx_max_no_tracks(SC_nuwro)){
            if (is_in_fv && SC_nuwro.npi0>0 && SC_nuwro.ccnc == 1) 
                h_pi0_e_nuwro_pass->Fill(SC_nuwro.pi0_e, weight);
        }
    }

    // Get the TTree
    _util.GetFile(f_mc_tune1, "../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_tune1.root"); // Get the run 1 MC file
    _util.GetTree(f_mc_tune1, mc_tree_tune1, "nuselection/NeutrinoSelectionFilter");

    SliceContainer SC_tune1;
    SC_tune1.Initialise(mc_tree_tune1, _util.k_mc, _util);

    int mc_tree_total_entries_tune1 = mc_tree_tune1->GetEntries();
    std::cout << "Total MC Events:         " << mc_tree_total_entries_tune1 << std::endl;

    // Event loop
    for (int ievent = 0; ievent < mc_tree_total_entries_tune1; ievent++){

        // See if we want to process all the events
        if (_util.num_events > 0){
            if (ievent >= _util.num_events) break;
        }

        // Alert the user
        if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
    
        // Get the entry in the tree
        mc_tree_tune1->GetEntry(ievent); 

        // Classify the event -- sets variable in the slice contianer
        SC_tune1.SliceClassifier(_util.k_mc);      // Classification of the event

        // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
        SC_tune1.SetThresholdEvent(_util.k_mc);

        // If the backtracked pdg of the leading shower is not an electron then alter classification
        SC_tune1.SetNonLdgShrEvent(_util.k_mc);
        
        SC_tune1.SliceInteractionType(_util.k_mc); // Genie interaction type
        SC_tune1.ParticleClassifier(_util.k_mc);   // The truth matched particle type of the leading shower
        SC_tune1.Pi0Classifier(_util.k_mc); 

        // Set derived variables in the slice container
        SC_tune1.SetSignal();                // Set the event as either signal or other
        SC_tune1.SetFakeData();              // Set the classifcation as data if fake data mode
        SC_tune1.SetTrueElectronThetaPhi();  // Set the true electron theta and phi variables
        SC_tune1.SetNuMIAngularVariables();  // Set the NuMI angular variables
        SC_tune1.CalibrateShowerEnergy();    // Divide the shower energy by 0.83 so it is done in one place

        // SC_tune1.SetPPFXCVWeight();

        bool is_in_fv = _util.in_fv(SC_tune1.true_nu_vtx_sce_x, SC_tune1.true_nu_vtx_sce_y, SC_tune1.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 

        double weight = _util.GetCVWeight(_util.k_mc, SC_tune1.weightSplineTimesTune, SC_tune1.ppfx_cv, SC_tune1.nu_e, SC_tune1.nu_pdg, is_in_fv, SC_tune1.interaction, SC_tune1.elec_e);
        
        // Pi0 Energy
        if (is_in_fv && SC_tune1.npi0>0 && SC_tune1.ccnc == 1) 
            h_pi0_e_tune1->Fill(SC_tune1.pi0_e, weight);


         // Require Passed cuts
        if (_scuts.swtrig(SC_tune1, _util.k_mc) &&
        _scuts.slice_id(SC_tune1)            &&
        _scuts.e_candidate(SC_tune1)         &&
        _scuts.in_fv(SC_tune1)               &&
        _scuts.contained_frac(SC_tune1)      &&
        _scuts.topo_score(SC_tune1)          &&
        _scuts.shr_cosmic_IP(SC_tune1)       &&
        _scuts.shower_score(SC_tune1)        &&
        _scuts.shr_hitratio(SC_tune1)        &&
        _scuts.shr_moliere_avg(SC_tune1)     &&
        _scuts.shr_dist_dEdx_max(SC_tune1)   &&
        _scuts.dEdx_max_no_tracks(SC_tune1)){
            if (is_in_fv && SC_tune1.npi0>0 && SC_tune1.ccnc == 1) 
                h_pi0_e_tune1_pass->Fill(SC_tune1.pi0_e, weight);
        }
    }





    // Now save the Pi0 Momentum plot
    _util.CreateDirectory("pi0");
    TCanvas * c = new TCanvas("c", "c", 500, 500);
    c->SetTopMargin(0.11);

    h_pi0_e_genie->SetStats(kFALSE);
    h_pi0_e_nuwro->SetStats(kFALSE);
    h_pi0_e_tune1->SetStats(kFALSE);

    _util.IncreaseLabelSize(h_pi0_e_genie, c);

    h_pi0_e_genie->SetLineColor(kAzure - 6);
    h_pi0_e_genie->SetLineWidth(2);
    h_pi0_e_genie->Draw("hist,E");

    h_pi0_e_nuwro->Scale(3.6248214);
    h_pi0_e_nuwro->SetLineColor(kRed+2);
    h_pi0_e_nuwro->SetLineWidth(2);
    h_pi0_e_nuwro->Draw("hist,E,same");

    h_pi0_e_tune1->Scale(3.8881733);
    h_pi0_e_tune1->SetLineColor(kGreen+2);
    h_pi0_e_tune1->SetLineWidth(2);
    h_pi0_e_tune1->Draw("hist,E,same");

    double yscale = h_pi0_e_genie->GetMaximum();
    if (h_pi0_e_nuwro->GetMaximum() > yscale)
        yscale = h_pi0_e_nuwro->GetMaximum();

    h_pi0_e_genie->SetMaximum(yscale*1.3);

    TLegend *leg = new TLegend(0.45, 0.89, 0.85, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->SetHeader(Form("GENIE/NuWro: %2.3f, GENIE/Tune 1: %2.3f", h_pi0_e_genie->Integral() / h_pi0_e_nuwro->Integral(), h_pi0_e_genie->Integral() / h_pi0_e_tune1->Integral()),"C");
    leg->AddEntry(h_pi0_e_genie, "GENIE v3", "l");
    leg->AddEntry(h_pi0_e_nuwro, "NuWro", "l");
    leg->AddEntry(h_pi0_e_tune1, "Tune 1", "l");
    leg->Draw();
    
    c->Print(Form("plots/run%s/pi0/h_pi0_comparison.pdf", _util.run_period));

    delete c;

    c = new TCanvas("c", "c", 500, 500);
    c->SetTopMargin(0.11);

    h_pi0_e_genie_pass->SetStats(kFALSE);
    h_pi0_e_nuwro_pass->SetStats(kFALSE);
    h_pi0_e_tune1_pass->SetStats(kFALSE);

    _util.IncreaseLabelSize(h_pi0_e_genie_pass, c);

    h_pi0_e_genie_pass->SetLineColor(kAzure - 6);
    h_pi0_e_genie_pass->SetLineWidth(2);
    h_pi0_e_genie_pass->Draw("hist,E");

    h_pi0_e_nuwro_pass->Scale(3.6248214);
    h_pi0_e_nuwro_pass->SetLineColor(kRed+2);
    h_pi0_e_nuwro_pass->SetLineWidth(2);
    h_pi0_e_nuwro_pass->Draw("hist,E,same");

    h_pi0_e_tune1_pass->Scale(3.8881733);
    h_pi0_e_tune1_pass->SetLineColor(kGreen+2);
    h_pi0_e_tune1_pass->SetLineWidth(2);
    h_pi0_e_tune1_pass->Draw("hist,E,same");

    yscale = h_pi0_e_genie_pass->GetMaximum();
    if (h_pi0_e_nuwro_pass->GetMaximum() > yscale)
        yscale = h_pi0_e_nuwro_pass->GetMaximum();

    h_pi0_e_genie_pass->SetMaximum(yscale*1.5);
    h_pi0_e_genie_pass->GetXaxis()->SetRangeUser(0, 1.5);

    TLegend *leg2 = new TLegend(0.45, 0.89, 0.85, 0.7);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg2->SetHeader(Form("GENIE/NuWro: %2.3f, GENIE/Tune 1: %2.3f", h_pi0_e_genie_pass->Integral() / h_pi0_e_nuwro_pass->Integral(), h_pi0_e_genie_pass->Integral() / h_pi0_e_tune1_pass->Integral()),"C");
    leg2->AddEntry(h_pi0_e_genie_pass, "GENIE v3", "l");
    leg2->AddEntry(h_pi0_e_nuwro_pass, "NuWro", "l");
    leg2->AddEntry(h_pi0_e_tune1_pass, "Tune 1", "l");
    leg2->Draw();
    
    c->Print(Form("plots/run%s/pi0/h_pi0_comparison_pass.pdf", _util.run_period));
 
    delete c;

}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareSelectedPi0(){

    gStyle->SetOptStat(0);

    // enums for legend
    enum enum_gens {k_cv, k_geniev3, k_tune1, k_nuwro, k_gen_MAX};
    std::vector<int> cols = {kBlack, kBlue+2, kOrange-2, kPink+1};
    
    std::vector<std::string> gen_names = {"CV", "geniev3", "tune1", "nuwro"};
    std::vector<double> sf = {1.0, 1.0, 3.8881733, 3.6248214};

    std::vector<TH1D*> h_pimass(k_gen_MAX);
    std::vector<TH1D*> h_shr_E(k_gen_MAX);
    std::vector<TH1D*> h_shr_cosbeta(k_gen_MAX);

    std::vector<TFile*> f(k_gen_MAX);

    TLegend *leg = new TLegend(0.45, 0.89, 0.85, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    std::string mode = "nc_pi0";
    // std::string mode = "numu_cc_pi0";

    for (unsigned int m = 0; m < gen_names.size(); m++){

        if (m == k_tune1)
            continue;
        
        if (m == k_cv)
            f.at(m) = TFile::Open(Form("files/nuexsec_mc_run1%s.root", ""), "READ");
        else
            f.at(m) = TFile::Open(Form("files/nuexsec_mc_run1_%s.root", gen_names.at(m).c_str()), "READ");

        if (m == k_cv)
            _util.GetHist(f.at(m), h_pimass.at(m), Form("pizero/h_pi0_mass_norm_%s",mode.c_str()));
        else
            _util.GetHist(f.at(m), h_pimass.at(m), Form("pizero/h_pi0_mass_%s", mode.c_str()));
        
        h_pimass.at(m)->Scale(sf.at(m));
        h_pimass.at(m)->SetLineColor(cols.at(m));
        h_pimass.at(m)->SetLineWidth(2);
        h_pimass.at(m)->SetTitle(Form("%s;#pi^{0} Mass [MeV];Entries", mode.c_str()));
        leg->AddEntry(h_pimass.at(m), gen_names.at(m).c_str(), "l");

        _util.GetHist(f.at(m), h_shr_E.at(m), Form("Stack/dEdx_max_no_tracks/%s/h_reco_shower_energy_cali_rebin_dEdx_max_no_tracks_%s", mode.c_str(), mode.c_str()));
        h_shr_E.at(m)->Scale(sf.at(m));
        h_shr_E.at(m)->SetLineColor(cols.at(m));
        h_shr_E.at(m)->SetLineWidth(2);
        h_shr_E.at(m)->Scale(1.0, "width");
        h_shr_E.at(m)->SetTitle(Form("%s;Shower Energy [GeV];Entries", mode.c_str()));

        _util.GetHist(f.at(m), h_shr_cosbeta.at(m), Form("Stack/dEdx_max_no_tracks/%s/h_reco_effective_cosangle_dEdx_max_no_tracks_%s", mode.c_str(), mode.c_str()));
        h_shr_cosbeta.at(m)->Scale(sf.at(m));
        h_shr_cosbeta.at(m)->SetLineColor(cols.at(m));
        h_shr_cosbeta.at(m)->SetLineWidth(2);
        h_shr_cosbeta.at(m)->Scale(1.0, "width");
        h_shr_cosbeta.at(m)->SetTitle(Form("%s;Shower cos#beta;Entries", mode.c_str()));


    }

    // Now Draw the histograms
    TCanvas * c, *c2, *c3; 
    
    for (unsigned int m = 0; m < gen_names.size(); m++){

        if (m == k_tune1)
            continue;

        
        if (m == k_cv){
            c = new TCanvas("c", "c", 500, 500);
            c->SetTopMargin(0.11);
            c2 = new TCanvas("c2", "c", 500, 500);
            c2->SetTopMargin(0.11);
            c3 = new TCanvas("c3", "c", 500, 500);
            c3->SetTopMargin(0.11);
        }

        c->cd();
        _util.IncreaseLabelSize(h_pimass.at(m), c);

        h_pimass.at(m)->Draw("hist,same");

        std::cout << gen_names.at(m) <<" Integral: "<< h_pimass.at(m)->Integral() << std::endl;

        double yscale = h_pimass.at(k_cv)->GetMaximum();
        h_pimass.at(k_cv)->SetMaximum(yscale*1.2);
        leg->Draw();

        c2->cd();
        _util.IncreaseLabelSize(h_shr_E.at(m), c2);

        h_shr_E.at(m)->Draw("hist,same");

        std::cout << gen_names.at(m) <<" Integral: "<< h_shr_E.at(m)->Integral() << std::endl;

        yscale = h_shr_E.at(k_cv)->GetMaximum();
        h_shr_E.at(k_cv)->SetMaximum(yscale*1.3);
        leg->Draw();

        c3->cd();
        _util.IncreaseLabelSize(h_shr_cosbeta.at(m), c3);

        h_shr_cosbeta.at(m)->Draw("hist,same");

        std::cout << gen_names.at(m) <<" Integral: "<< h_shr_cosbeta.at(m)->Integral() << std::endl;

        yscale = h_shr_cosbeta.at(k_cv)->GetMaximum();
        h_shr_cosbeta.at(k_cv)->SetMaximum(yscale*1.3);
        leg->Draw();

        if (m == k_gen_MAX-1){

            c->Print(Form("plots/run%s/pi0/h_%s_mass_comparison.pdf", _util.run_period, mode.c_str()));
            c2->Print(Form("plots/run%s/pi0/h_%s_energy_comparison.pdf", _util.run_period, mode.c_str()));
            c3->Print(Form("plots/run%s/pi0/h_%s_cosbeta_comparison.pdf", _util.run_period, mode.c_str()));
            delete c;
            delete c2;
        }
    }

    

}
// -----------------------------------------------------------------------------
void UtilityPlotter::Compare1DFluxGeneratedEvents(){

    gStyle->SetOptStat(0);

    std::string flav = "nuebar";
    std::string var = "elec_e"; // nu_e, elec_e, cosbeta

    std::string title = "#nu_{e}";
    if (flav == "nuebar")
        title = "#bar{#nu}_{e}";

    // Load in the Nue file
    TFile *f_nue;    
    TTree *t_nue;
    
    
    if (flav == "nue"){
        f_nue = TFile::Open("../ntuples/genie_gen_nue.gst.root", "READ");
        _util.GetTree(f_nue, t_nue, "gst");
    }
    else {
       f_nue = TFile::Open("../ntuples/genie_gen_nuebar.gst.root", "READ");
       _util.GetTree(f_nue, t_nue, "gst");
    }

    TFile *f_xsec   = TFile::Open("files/trees/nuexsec_selected_tree_mc_run1.root");
    TTree *t_xsec;
    _util.GetTree(f_xsec, t_xsec, "mc_nue_tree");
    
    std::string signal_def = "cc==1 && Ev>0.06 && El>0.12";
    std::string signal_def2;
    
    if (flav == "nue")
        signal_def2 = "ppfx_cv*(nu_pdg == 12)";
    else
        signal_def2 = "ppfx_cv*(nu_pdg == -12)";

    std::string query;
    
    if (var == "cosbeta")
        query = "(pxv*pxl + pyv*pyl + pzv*pzl) / (sqrt(pxv*pxv + pyv*pyv + pzv*pzv) * sqrt(pxl*pxl + pyl*pyl + pzl*pzl) )";
    else if (var == "elec_e")
        query = "El";
    else
        query = "Ev";
        

    TCanvas * c = new TCanvas("c", "c", 500, 500);
    TH1D *htemp, *htemp2;
    
    if (var == "cosbeta"){
        htemp  = new TH1D("htemp", Form("%s;cos#beta;", title.c_str()), 100, -1, 1);
        htemp2 = new TH1D("htemp2",Form("%s;cos#beta;", title.c_str()), 100, -1, 1);
    }
    else if (var == "elec_e"){
        htemp   = new TH1D("htemp", Form("%s;Electron Energy [GeV];", title.c_str()), 100, 0, 5);
        htemp2  = new TH1D("htemp2",Form("%s;Electron Energy [GeV];", title.c_str()), 100, 0, 5);
    }
    else {
        htemp   = new TH1D("htemp", Form("%s;Neutrino Energy [GeV];", title.c_str()), 100, 0, 5);
        htemp2  = new TH1D("htemp2",Form("%s;Neutrino Energy [GeV];", title.c_str()), 100, 0, 5);
    }

    // Draw the Query -- adjust by query type
    t_nue->Draw(Form("%s >> htemp", query.c_str()), signal_def.c_str());
    
    if (var == "cosbeta")
        t_xsec->Draw(Form("%s >> htemp2", "cos_true_effective_angle"), signal_def2.c_str());
    else if (var == "elec_e")
        t_xsec->Draw(Form("%s >> htemp2", "elec_e"), signal_def2.c_str());
    else
        t_xsec->Draw(Form("%s >> htemp2", "true_energy"), signal_def2.c_str());

    TLegend *leg = new TLegend(0.45, 0.89, 0.85, 0.7);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(htemp2, "MC", "l");
    leg->AddEntry(htemp, "Gevgen 1D Flux", "l");


    // Draw the histogram
    htemp->SetLineWidth(2);
    htemp->SetLineColor(kBlack);
    htemp2->SetLineWidth(2);
    htemp2->SetLineColor(kBlue+2);
    
    htemp->Scale(htemp2->Integral()/htemp->Integral());

    _util.IncreaseLabelSize(htemp, c);
    
    htemp->Draw("hist,E");
    htemp2->Draw("hist,E,same");

    leg->Draw();

    TPaveText *pt;

    pt = new TPaveText(0.3215, 0.936, 0.3215, 0.936, "NDC");
    pt->AddText("Area Normalised");
    pt->SetTextColor(kGreen + 2);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();

    c->Print(Form("plots/run%s/Truth/1D_flux_comparisons_%s_%s.pdf", _util.run_period, flav.c_str(),var.c_str() ));
}
// -----------------------------------------------------------------------------
void UtilityPlotter::CompareXsecPi0Tunings(){


   gStyle->SetOptStat(0);

    // Create a vector for the models
    std::vector<std::string> models = {
        "Input",
        "nopi0tune"
    };

    // enums for the models
    enum enum_models {
        k_model_input,
        k_model_nopi0tune,
        k_MODEL_MAX
    };

    // Load in the cross section output
    TFile *fxsec = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "READ");

    TH2D* h_temp_2D;
    TH1D* h_temp;

    // Total Covariance Matrix
    h_temp_2D = (TH2D*)fxsec->Get(Form("%s/wiener/h_cov_tot_dataxsec_reco",_util.xsec_var));
    TH2D* h_cov_reco = (TH2D*)h_temp_2D->Clone();
    h_cov_reco->SetDirectory(0);

    fxsec->Close();

    // Now Get the Models
    // Load in the cross section output
    fxsec = TFile::Open(Form("files/crosssec_run%s.root ", _util.run_period), "READ");

    std::vector<TH1D*> h_true_model(k_MODEL_MAX);
    std::vector<TH1D*> h_data_model(k_MODEL_MAX);
    std::vector<TH2D*> h_cov_diag(k_MODEL_MAX);
    std::vector<TH2D*> h_response_model(k_MODEL_MAX);
    std::vector<TH1D*> h_unf_model(k_MODEL_MAX);

    // Make the plot
    TCanvas *c = new TCanvas("c", "c", 500, 500);
    
    // Create the legend
    TLegend *leg = new TLegend(0.35, 0.50, 0.70, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // Loop over each model
    for (unsigned int m = 0; m < models.size(); m++){

        // Get true model xsec
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run%s_CV_0_%s_mc_xsec", models.at(m).c_str(),_util.vars.at(k_var_trueX).c_str(), _util.run_period, _util.vars.at(k_var_trueX).c_str()));
        h_true_model.at(m) = (TH1D*)h_temp->Clone();

        // Get data cross section extracted for model 
        h_temp  = (TH1D*)fxsec->Get(Form("%s/%s/h_run%s_CV_0_%s_data_xsec", models.at(m).c_str(), _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
        h_data_model.at(m) = (TH1D*)h_temp->Clone();

        // Get the response matrix for model
        h_temp_2D  = (TH2D*)fxsec->Get(Form("%s/%s/h_run%s_%s_0_smearing", models.at(m).c_str(), _util.vars.at(k_var_trueX).c_str(), _util.run_period, models.at(m).c_str()));
        h_response_model.at(m) = (TH2D*)h_temp_2D->Clone();

        // Flip the response matrix
        // Loop over rows
        for (int row=0; row<h_temp_2D->GetXaxis()->GetNbins()+2; row++) {

            for (int col=0; col<h_temp_2D->GetYaxis()->GetNbins()+2; col++){
                h_response_model.at(m)->SetBinContent(col, row, h_temp_2D->GetBinContent(row, col));          
            }
        }

        
        h_cov_diag.at(m) = (TH2D*)h_cov_reco->Clone();
        
        // Initialise the WienerSVD class
        WienerSVD _wSVD;
        _wSVD.Initialise(_util);
        _wSVD.DoUnfolding(2, 0, h_true_model.at(m), h_data_model.at(m), h_response_model.at(m), h_cov_diag.at(m));

        h_unf_model.at(m) = (TH1D*)_wSVD.unf->Clone(Form("test_%s", models.at(m).c_str()));

        for (int bin = 1; bin < h_unf_model.at(m)->GetNbinsX()+1; bin++){
            double err = _wSVD.unfcov->GetBinContent(bin, bin);
            h_unf_model.at(m)->SetBinError(bin, std::sqrt(err));
        }
        
        // Set the line colours
        if (m == k_model_input){
            h_unf_model.at(m)->SetLineColor(kBlack);
            h_unf_model.at(m)->SetMarkerStyle(20);
            h_unf_model.at(m)->SetMarkerSize(0.5);
        }    
        

        if (m == k_model_nopi0tune)    h_unf_model.at(m)->SetLineColor(kPink+1);

        
        h_unf_model.at(m)->Scale(1.0, "width");
        
        if (m == k_model_input){
            _util.IncreaseLabelSize( h_unf_model.at(m), c);
            gPad->SetLeftMargin(0.20);
            c->SetBottomMargin(0.15);
            h_unf_model.at(m)->SetTitle(_util.var_labels_xsec.at(k_var_trueX).c_str());
            h_unf_model.at(m)->Draw("E,same");
            h_unf_model.at(m)->GetYaxis()->SetTitleOffset(1.4);
        }
        else {
            h_unf_model.at(m)->Draw("hist,same");
        }

        h_unf_model.at(m)->SetLineWidth(2);

        h_unf_model.at(m)->SetMinimum(0);

        
        if (m == k_model_input) {
            models.at(m) = "CV";
            leg->AddEntry(h_unf_model.at(m), "Data (Stat. + Sys.)", "ep");
        }
        if (m == k_model_nopi0tune)      leg->AddEntry(h_unf_model.at(m), "Data with #pi^{0} Tune", "l");
        

        delete _wSVD.smear;
        delete _wSVD.wiener;
        delete _wSVD.unfcov;
        delete _wSVD.unf;
        delete _wSVD.diff;
        delete _wSVD.bias;
        delete _wSVD.bias2;
        delete _wSVD.fracError;
        delete _wSVD.absError;
        delete _wSVD.MSE;
        delete _wSVD.MSE2;
    }

    h_unf_model.at(k_model_input)->Draw("E,same");
    gStyle->SetLegendTextSize(0.04);
    leg->Draw();

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.52, 0.92, 0.52, 0.92);

    c->Print(Form("plots/run%s/Models/%s/run%s_UnfoldedDataComparison_%s_pi0tune.pdf", _util.run_period, _util.xsec_var, _util.run_period, _util.xsec_var));
    delete c;




}
// -----------------------------------------------------------------------------
void UtilityPlotter::CalcFluxCovarianceHP(){

    // enums for the models
    enum enum_flav {
        k_numu,
        k_numubar,
        k_nue,
        k_nuebar,
        k_FLAV_MAX
    };

    std::string mode = "2D";
    mode = "1D"; // Only study the histograms in terms of energy


    std::vector<TH1D*> hist_unwrap_stitch;
    TH1D *hist_unwrap_stitch_CV;
    std::vector<TH1D*> hist_unwrap_CV(k_FLAV_MAX);

    // Call function here
    GetStitchedUniverses("HP", mode, hist_unwrap_stitch, hist_unwrap_stitch_CV, hist_unwrap_CV, 1 );
    
    // Draw vertical lines to help the eye
    std::vector<TLine*> line(6);

    int numu_bin    = hist_unwrap_CV.at(k_numu)->GetNbinsX() + 1;
    int numubar_bin = hist_unwrap_CV.at(k_numu)->GetNbinsX() + hist_unwrap_CV.at(k_numubar)->GetNbinsX() + 1;
    int nue_bin     = hist_unwrap_CV.at(k_numu)->GetNbinsX() + hist_unwrap_CV.at(k_numubar)->GetNbinsX() + hist_unwrap_CV.at(k_nue)->GetNbinsX() + 1;
    int nuebar_bin  = hist_unwrap_CV.at(k_numu)->GetNbinsX() + hist_unwrap_CV.at(k_numubar)->GetNbinsX() + hist_unwrap_CV.at(k_nue)->GetNbinsX() + hist_unwrap_CV.at(k_nuebar)->GetNbinsX() + 1;

    line.at(0) = new TLine(numu_bin,    1, numu_bin,    nuebar_bin);
    line.at(1) = new TLine(numubar_bin, 1, numubar_bin, nuebar_bin);
    line.at(2) = new TLine(nue_bin,     1, nue_bin,     nuebar_bin);
    line.at(3) = new TLine(1,numu_bin,    nuebar_bin, numu_bin);
    line.at(4) = new TLine(1,numubar_bin, nuebar_bin, numubar_bin);
    line.at(5) = new TLine(1,nue_bin,     nuebar_bin, nue_bin);


    // Create the covariance matrix
    int n_bins = hist_unwrap_stitch.at(0)->GetNbinsX();
    TH2D* h_cov = new TH2D("", "Covariance Matrix ;Bin i; Bin j",  n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);

    _util.CalcCovariance(hist_unwrap_stitch, hist_unwrap_stitch_CV, h_cov);

    gStyle->SetOptStat(0);
    
    // gStyle->SetPalette(kBlueGreenYellow);

    TH2D* h_cor = (TH2D*)h_cov->Clone();
    TH2D* h_frac_cov = (TH2D*)h_cov->Clone();
    _util.CalcCorrelation(hist_unwrap_stitch_CV, h_cov, h_cor);
    _util.CalcCFracCovariance(hist_unwrap_stitch_CV, h_frac_cov);

    TCanvas *c = new TCanvas("", "", 500, 500);
    h_cov->Draw("colz");
    for (unsigned int i = 0; i< line.size(); i++){
        line.at(i)->SetLineColor(kRed+2);
        line.at(i)->SetLineWidth(4);
        line.at(i)->Draw();
    }
    c->Print("covariance.pdf");

    TCanvas *c2 = new TCanvas("", "", 700, 700);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.14);
    gPad->SetTopMargin(0.14);
    gPad->SetBottomMargin(0.14);
    h_cor->SetTitle("Correlation Matrix ;Bin i; Bin j");
    h_cor->Draw("colz");
    for (unsigned int i = 0; i< line.size(); i++){
        line.at(i)->SetLineColor(kRed+2);
        line.at(i)->SetLineWidth(4);
        line.at(i)->Draw();
    }
    c2->Print("correlation.pdf");

    TCanvas *c3 = new TCanvas("", "", 500, 500);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.14);
    gPad->SetTopMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // h_cor->SetMaximum(0.3);
    // h_cor->SetMinimum(-0.3);
    h_frac_cov->SetTitle("Fractional Covariance Matrix ;Bin i; Bin j");
    h_frac_cov->Draw("colz");
    // gPad->SetLogz();
    for (unsigned int i = 0; i< line.size(); i++){
        line.at(i)->SetLineColor(kRed+2);
        line.at(i)->SetLineWidth(4);
        line.at(i)->Draw();
    }
    c3->Print("fraction_cov.pdf");



   


}
// -----------------------------------------------------------------------------
void UtilityPlotter::CalcFluxCovarianceBeamline(){

    // enums for the models
    enum enum_flav {
        k_numu,
        k_numubar,
        k_nue,
        k_nuebar,
        k_FLAV_MAX
    };

    std::string mode = "2D";
    // mode = "1D"; // Only study the histograms in terms of energy

    // std::vector<int> index = {1,3,5,7,9,11,13,15,17,19};
    std::vector<int> index = {9};

    std::vector<std::vector<TH1D*>> hist_unwrap_stitch;
    hist_unwrap_stitch.resize(index.size());

    std::vector<TH1D*> hist_unwrap_stitch_CV;
    hist_unwrap_stitch_CV.resize(index.size());

    std::vector<std::vector<TH1D*>> hist_unwrap_CV(k_FLAV_MAX);
    hist_unwrap_CV.resize(index.size());
    for (unsigned int i = 0; i < hist_unwrap_CV.size(); i++)
        hist_unwrap_CV.at(i).resize(k_FLAV_MAX);

    for (unsigned int i = 0; i < index.size(); i++){
        GetStitchedUniverses("beamline", mode, hist_unwrap_stitch.at(i), hist_unwrap_stitch_CV.at(i), hist_unwrap_CV.at(i), index.at(i) );
    }

    // Draw vertical lines to help the eye
    std::vector<TLine*> line(6);

    int numu_bin    = hist_unwrap_CV.at(0).at(k_numu)->GetNbinsX() + 1;
    int numubar_bin = hist_unwrap_CV.at(0).at(k_numu)->GetNbinsX() + hist_unwrap_CV.at(0).at(k_numubar)->GetNbinsX() + 1;
    int nue_bin     = hist_unwrap_CV.at(0).at(k_numu)->GetNbinsX() + hist_unwrap_CV.at(0).at(k_numubar)->GetNbinsX() + hist_unwrap_CV.at(0).at(k_nue)->GetNbinsX() + 1;
    int nuebar_bin  = hist_unwrap_CV.at(0).at(k_numu)->GetNbinsX() + hist_unwrap_CV.at(0).at(k_numubar)->GetNbinsX() + hist_unwrap_CV.at(0).at(k_nue)->GetNbinsX() + hist_unwrap_CV.at(0).at(k_nuebar)->GetNbinsX() + 1;

    line.at(0) = new TLine(numu_bin,    1, numu_bin,    nuebar_bin);
    line.at(1) = new TLine(numubar_bin, 1, numubar_bin, nuebar_bin);
    line.at(2) = new TLine(nue_bin,     1, nue_bin,     nuebar_bin);
    line.at(3) = new TLine(1,numu_bin,    nuebar_bin, numu_bin);
    line.at(4) = new TLine(1,numubar_bin, nuebar_bin, numubar_bin);
    line.at(5) = new TLine(1,nue_bin,     nuebar_bin, nue_bin);


    // Create the covariance matrix
    int n_bins = hist_unwrap_stitch.at(0).at(0)->GetNbinsX();
    
    std::vector<TH2D*> h_cov_v(index.size());
    
    for (unsigned int i = 0; i < index.size(); i++){
        h_cov_v.at(i) = new TH2D(Form("cov_%d", i), "Covariance Matrix ;Bin i; Bin j",  n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);

        _util.CalcCovariance(hist_unwrap_stitch.at(i), hist_unwrap_stitch_CV.at(i), h_cov_v.at(i));
    }

    TH2D* h_cov = (TH2D*)h_cov_v.at(0)->Clone();
    for (unsigned int i = 1; i < index.size(); i++){
        h_cov->Add(h_cov_v.at(i));
    }

    gStyle->SetOptStat(0);
    
    // gStyle->SetPalette(kBlueGreenYellow);

    TH2D* h_cor = (TH2D*)h_cov->Clone();
    TH2D* h_frac_cov = (TH2D*)h_cov->Clone();
    _util.CalcCorrelation(hist_unwrap_stitch_CV.at(0), h_cov, h_cor);
    _util.CalcCFracCovariance(hist_unwrap_stitch_CV.at(0), h_frac_cov);

    TCanvas *c = new TCanvas("", "", 500, 500);
    h_cov->Draw("colz");
    for (unsigned int i = 0; i< line.size(); i++){
        line.at(i)->SetLineColor(kRed+2);
        line.at(i)->SetLineWidth(4);
        line.at(i)->Draw();
    }
    c->Print("covariance.pdf");

    TCanvas *c2 = new TCanvas("", "", 700, 700);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.14);
    gPad->SetTopMargin(0.14);
    gPad->SetBottomMargin(0.14);
    h_cor->SetTitle("Correlation Matrix ;Bin i; Bin j");
    h_cor->Draw("colz");
    for (unsigned int i = 0; i< line.size(); i++){
        line.at(i)->SetLineColor(kRed+2);
        line.at(i)->SetLineWidth(4);
        line.at(i)->Draw();
    }
    c2->Print("correlation.pdf");

    TCanvas *c3 = new TCanvas("", "", 500, 500);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.14);
    gPad->SetTopMargin(0.14);
    gPad->SetBottomMargin(0.14);
    // h_cor->SetMaximum(0.3);
    // h_cor->SetMinimum(-0.3);
    h_frac_cov->SetTitle("Fractional Covariance Matrix ;Bin i; Bin j");
    h_frac_cov->Draw("colz");
    // gPad->SetLogz();
    for (unsigned int i = 0; i< line.size(); i++){
        line.at(i)->SetLineColor(kRed+2);
        line.at(i)->SetLineWidth(4);
        line.at(i)->Draw();
    }
    c3->Print("fraction_cov.pdf");

}
// -----------------------------------------------------------------------------
void UtilityPlotter::GetStitchedUniverses(std::string constraint, std::string mode, std::vector<TH1D*> &hist_unwrap_stitch, TH1D* &hist_unwrap_stitch_CV, std::vector<TH1D*> &hist_unwrap_CV, int index){

    // Create a vector for the models
    std::vector<std::string> flav = {
        "numu",
        "numubar",
        "nue",
        "nuebar"
    };

    // enums for the models
    enum enum_flav {
        k_numu,
        k_numubar,
        k_nue,
        k_nuebar,
        k_FLAV_MAX
    };

    int nuniverses = 600;
    if (constraint =="beamline")
        nuniverses = 2;

    // Load in the flux histogram file
    TFile * f = TFile::Open("Systematics/output_fhc_uboone_run0.root", "READ");

    TFile * f2, *f3;
    if (constraint == "beamline"){
        f2 = TFile::Open(Form("Systematics//beamline/FHC/output_uboone_fhc_run%d.root", index), "READ");
        f3 = TFile::Open(Form("Systematics//beamline/FHC/output_uboone_fhc_run%d.root", index+1), "READ");
    }


    // Create the 2D histograms to get from file
    std::vector<std::vector<TH2D*>> hist;
    hist.resize(k_FLAV_MAX);
    for (unsigned int h = 0; h < hist.size(); h++ ){
        hist.at(h).resize(nuniverses);
    }

    // Create the vector of histograms for unwrapping
    std::vector<std::vector<TH1D*>> hist_unwrap;
    hist_unwrap.resize(k_FLAV_MAX);
    for (unsigned int h = 0; h < hist_unwrap.size(); h++ ){
        hist_unwrap.at(h).resize(nuniverses);
    }

    std::vector<TH2D*> hist_2D_CV(k_FLAV_MAX);

    // 2D histogram case
    if (mode == "2D"){
        // Now get the histograms
        std::cout << "Reading in the histograms from file..." << std::endl;
        for (unsigned int h = 0; h < hist.size(); h++ ){
            
            // Get the histograms
            for (unsigned int u = 0; u < hist.at(h).size(); u++ ){
                
                if (constraint == "HP")
                    hist.at(h).at(u) = (TH2D*)f->Get(Form("%s/Multisims/%s_ppfx_ms_UBPPFX_Uni_%d_AV_TPC_2D", flav.at(h).c_str(), flav.at(h).c_str(), u));

                if (constraint == "beamline"){
                    hist.at(h).at(0) = (TH2D*)f2->Get(Form("%s/Detsmear/%s_CV_AV_TPC_2D", flav.at(h).c_str(), flav.at(h).c_str()));
                    hist.at(h).at(1) = (TH2D*)f3->Get(Form("%s/Detsmear/%s_CV_AV_TPC_2D", flav.at(h).c_str(), flav.at(h).c_str()));
                }
            }
            
            // Get the CV
            TH2D* h_temp;
            h_temp= (TH2D*)f->Get(Form("%s/Detsmear/%s_CV_AV_TPC_2D", flav.at(h).c_str(), flav.at(h).c_str()));
            hist_2D_CV.at(h) = (TH2D*)h_temp->Clone();


        }
    }

    std::cout << "Successfully read in the histograms from file" << std::endl;

    // --
    // Unwrap the histograms
    

    // Loop over the flavours
    for (unsigned int h = 0; h < hist.size(); h++ ){
        
        // Loop over the universes
        for (unsigned int u = 0; u < hist.at(h).size(); u++ ){

            if (mode == "2D"){
                // Create the unwrapped histograms
                const int nBinsEnu = hist.at(h).at(0)->GetXaxis()->GetNbins(); // Enu
                const int nBinsTh  = hist.at(h).at(0)->GetYaxis()->GetNbins(); // Theta
                hist_unwrap.at(h).at(u) = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh);
                
                // Fill the unwrapped histograms
                int counter{0};
                for (int i=1; i<hist.at(h).at(u)->GetXaxis()->GetNbins()+1; i++) { // Loop over rows
                    for (int j=1; j<hist.at(h).at(u)->GetYaxis()->GetNbins()+1; j++){// Loop over columns
                        counter++;
                        hist_unwrap.at(h).at(u)->SetBinContent(counter, hist.at(h).at(u)->GetBinContent(i , j)  );
                    }
                }
            }
            // 1D histogram case
            else {

                if (constraint == "HP"){
                    hist_unwrap.at(h).at(u) = (TH1D*)f->Get(Form("%s/Multisims/%s_ppfx_ms_UBPPFX_Uni_%d_AV_TPC", flav.at(h).c_str(), flav.at(h).c_str(), u));
                }

                if (constraint == "beamline"){
                    hist_unwrap.at(h).at(0) = (TH1D*)f2->Get(Form("%s/Detsmear/%s_CV_AV_TPC", flav.at(h).c_str(), flav.at(h).c_str()));
                    hist_unwrap.at(h).at(1) = (TH1D*)f3->Get(Form("%s/Detsmear/%s_CV_AV_TPC", flav.at(h).c_str(), flav.at(h).c_str()));
                }
            
            }

        }

        // CV
        if (mode == "2D"){
            // Create the unwrapped histograms
            const int nBinsEnu = hist.at(h).at(0)->GetXaxis()->GetNbins(); // Enu
            const int nBinsTh  = hist.at(h).at(0)->GetYaxis()->GetNbins(); // Theta
            hist_unwrap_CV.at(h) = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh);

            // Fill the CV unwrapped histograms
            int counter{0};
            for (int i=1; i<hist_2D_CV.at(h)->GetXaxis()->GetNbins()+1; i++) { // Loop over rows
                for (int j=1; j<hist_2D_CV.at(h)->GetYaxis()->GetNbins()+1; j++){ // Loop over columns
                    counter++;
                    hist_unwrap_CV.at(h)->SetBinContent(counter, hist_2D_CV.at(h)->GetBinContent(i , j)  );
                }
            }

        }
        else {
            
            TH1D* h_temp;
            h_temp= (TH1D*)f->Get(Form("%s/Detsmear/%s_CV_AV_TPC", flav.at(h).c_str(), flav.at(h).c_str()));
            hist_unwrap_CV.at(h) = (TH1D*)h_temp->Clone();
        }


    }

    std::cout << "Successfully unwrapped the histograms" << std::endl;

    // Now stitch together all the histograms
    
    hist_unwrap_stitch.resize(nuniverses);
    // Loop over universes
    for (int u = 0; u < nuniverses; u++ ){
        
        // Create the stitched histogram
        const int nBins = hist_unwrap.at(k_numu).at(0)->GetNbinsX() + hist_unwrap.at(k_numubar).at(0)->GetNbinsX() + hist_unwrap.at(k_nue).at(0)->GetNbinsX() + hist_unwrap.at(k_nuebar).at(0)->GetNbinsX();
        hist_unwrap_stitch.at(u) = new TH1D("", "",nBins, 0, nBins);

        std::vector<double> values;

        // Loop over histograms for each flavour and get the bin contents.
        for (unsigned int h = 0; h < hist.size(); h++ ){

            for (int bin=1; bin<hist_unwrap.at(h).at(u)->GetXaxis()->GetNbins()+1; bin++){
                values.push_back( hist_unwrap.at(h).at(u)->GetBinContent(bin) );
            }
        }

        // Set the bin content of the stitched histogram
        for (unsigned int bin=1; bin<values.size()+1; bin++) { // Loop over rows
    
            hist_unwrap_stitch.at(u)->SetBinContent(bin, values.at(bin-1));
        }

    }

    // Create the stitched histogram
    const int nBins = hist_unwrap_CV.at(k_numu)->GetNbinsX() + hist_unwrap_CV.at(k_numubar)->GetNbinsX() + hist_unwrap_CV.at(k_nue)->GetNbinsX() + hist_unwrap_CV.at(k_nuebar)->GetNbinsX();
    hist_unwrap_stitch_CV = new TH1D("", "",nBins, 0, nBins);

    std::vector<double> values;

    // Loop over histograms for each flavour and get the bin contents.
    for (unsigned int h = 0; h < hist.size(); h++ ){

        for (int bin=1; bin<hist_unwrap_CV.at(h)->GetXaxis()->GetNbins()+1; bin++){
            values.push_back( hist_unwrap_CV.at(h)->GetBinContent(bin) );
        }
    }

    // Set the bin content of the stitched histogram
    for (unsigned int bin=1; bin<values.size()+1; bin++) { // Loop over rows

        hist_unwrap_stitch_CV->SetBinContent(bin, values.at(bin-1));
    }

}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotParentEventRates(std::string type){


    TFile *f_standard;
    TTree *t_standard;
    
    TFile *f_intrinsic;
    TTree *t_intrinsic;

    std::string title;

    std::vector<std::string> par_pdg_str = {"#pi^{+}", "#pi^{-}", "#mu^{+}", "#mu^{-}","K^{+}","K^{-}", "K^{0}_{L}", "Total" };
    std::vector<std::string> par_pdg_str2 = {"pip", "pim", "mup", "mum","kp","km", "kl", "all" };
    enum par_enum {k_pi_plus = 211, k_pi_minus = -211, k_mu_plus = -13, k_mu_minus = 13, k_K_plus = 321, k_K_minus = -321, k_k0 = 130, k_all = 1000};
    std::vector<int> par_pdg = {211, -211, -13, 13, 321, -321, 130, 0};
    std::vector<int> cols = {42, kMagenta+2, 30, 38, 28, 36, 1001, kBlack };

    std::vector<std::string> flav = {"#nu_{#mu}", "#bar{#nu}_{#mu}", "#nu_{e}", "#bar{#nu}_{e}"};
    std::vector<std::string> flav_norm = {"numu", "numubar", "nue", "nuebar"};


    std::vector<std::string> q_numu(par_pdg.size());
    std::vector<std::string> q_numubar(par_pdg.size());
    std::vector<std::string> q_nue(par_pdg.size());
    std::vector<std::string> q_nuebar(par_pdg.size());
    std::string q_all;

    // Load in the TFiles
    if (type == "dk2nu"){
        f_standard = TFile::Open("../ntuples/neutrinoselection_filt_run1_overlay_newtune.root" , "READ");
        t_standard = (TTree*)f_standard->Get("nuselection/NeutrinoSelectionFilter");
        
        f_intrinsic = TFile::Open("../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root" , "READ");
        t_intrinsic = (TTree*)f_intrinsic->Get("nuselection/NeutrinoSelectionFilter");

        title = "dk2nu";
        for (unsigned int pdg = 0; pdg < par_pdg.size(); pdg++){
            q_numu.at(pdg)    = Form("ppfx_cv*(ccnc==0 && nu_pdg == 14  && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91 && nu_parent_pdg == %d)", par_pdg.at(pdg));
            q_numubar.at(pdg) = Form("ppfx_cv*(ccnc==0 && nu_pdg == -14 && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91 && nu_parent_pdg == %d)", par_pdg.at(pdg));
            q_nue.at(pdg)     = Form("ppfx_cv*(ccnc==0 && nu_pdg == 12  && nu_parent_pdg == %d)", par_pdg.at(pdg));
            q_nuebar.at(pdg)  = Form("ppfx_cv*(ccnc==0 && nu_pdg == -12 && nu_parent_pdg == %d)", par_pdg.at(pdg));
        
            if (par_pdg.at(pdg) == 0){
                q_numu.at(pdg)    = "ppfx_cv*(ccnc==0 && nu_pdg == 14  && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91)";
                q_numubar.at(pdg) = "ppfx_cv*(ccnc==0 && nu_pdg == -14 && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91)";
                q_nue.at(pdg)     = "ppfx_cv*(ccnc==0 && nu_pdg == 12 )";
                q_nuebar.at(pdg)  = "ppfx_cv*(ccnc==0 && nu_pdg == -12)";
            }
        }
    }
    else if (type == "flugg"){

        f_standard = TFile::Open("../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_FLUGG.root" , "READ");
        t_standard = (TTree*)f_standard->Get("nuselection/NeutrinoSelectionFilter");

        f_intrinsic = TFile::Open("../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_FLUGG_intrinsic.root" , "READ");
        t_intrinsic = (TTree*)f_intrinsic->Get("nuselection/NeutrinoSelectionFilter");

        title = "flugg";
        for (unsigned int pdg = 0; pdg < par_pdg.size(); pdg++){
            q_numu.at(pdg)    = Form("(ccnc==0 && nu_pdg == 14  && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91 && nu_parent_pdg == %d)", par_pdg.at(pdg));
            q_numubar.at(pdg) = Form("(ccnc==0 && nu_pdg == -14 && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91 && nu_parent_pdg == %d)", par_pdg.at(pdg));
            q_nue.at(pdg)     = Form("(ccnc==0 && nu_pdg == 12  && nu_parent_pdg == %d)", par_pdg.at(pdg));
            q_nuebar.at(pdg)  = Form("(ccnc==0 && nu_pdg == -12 && nu_parent_pdg == %d)", par_pdg.at(pdg));
        
            if (par_pdg.at(pdg) == 0){
                q_numu.at(pdg)    = "(ccnc==0 && nu_pdg == 14  && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91)";
                q_numubar.at(pdg) = "(ccnc==0 && nu_pdg == -14 && true_nu_vtx_x >= -1.55 && true_nu_vtx_x <= 254.8 && true_nu_vtx_y >= -115.53 && true_nu_vtx_y <= 117.47 && true_nu_vtx_z >= 0.11 && true_nu_vtx_z <= 1036.91)";
                q_nue.at(pdg)     = "(ccnc==0 && nu_pdg == 12 )";
                q_nuebar.at(pdg)  = "(ccnc==0 && nu_pdg == -12)";
            }
        }
    }
    else {
        return;
    }
    
    std::vector<TH1D*> h_numu(par_pdg.size()); 
    std::vector<TH1D*> h_numubar(par_pdg.size()); 
    std::vector<TH1D*> h_nue(par_pdg.size()); 
    std::vector<TH1D*> h_nuebar(par_pdg.size());
    
    for (unsigned int pdg = 0; pdg < par_pdg.size(); pdg++){
        h_numu.at(pdg)    = new TH1D(Form("h_numu_%s", par_pdg_str.at(pdg).c_str()),   Form("%s #nu_{#mu};Neutrino Energy [GeV];Entries",        title.c_str()),50, 0, 5);
        h_numubar.at(pdg) = new TH1D(Form("h_numubar_%s", par_pdg_str.at(pdg).c_str()),Form("%s #bar{#nu}_{#mu};Neutrino Energy [GeV];Entries", title.c_str()), 35, 0, 5);
        h_nue.at(pdg)     = new TH1D(Form("h_nue_%s", par_pdg_str.at(pdg).c_str()),    Form("%s #nu_{e};Neutrino Energy [GeV];Entries",         title.c_str()), 50, 0, 5);
        h_nuebar.at(pdg)  = new TH1D(Form("h_nuebar_%s", par_pdg_str.at(pdg).c_str()), Form("%s #bar{#nu}_{e};Electron Energy [GeV];Entries",   title.c_str()), 30, 0, 5);
    }
    
    std::vector<TLegend*> leg(flav.size());
    for (unsigned int f = 0; f < leg.size(); f++){
        leg.at(f) = new TLegend(0.64, 0.45, 0.84, 0.9);
        leg.at(f)->SetBorderSize(0);
        leg.at(f)->SetFillStyle(0);
        leg.at(f)->SetTextSize(0.05);
    }

    for (unsigned int pdg = 0; pdg < par_pdg.size(); pdg++){
        std::cout << "On Query: " << par_pdg.at(pdg) << std::endl;

        t_standard->Draw(Form("nu_e >> h_numu_%s",    par_pdg_str.at(pdg).c_str()), q_numu.at(pdg).c_str());
        t_standard->Draw(Form("nu_e >> h_numubar_%s", par_pdg_str.at(pdg).c_str()), q_numubar.at(pdg).c_str());
        t_intrinsic->Draw(Form("nu_e >> h_nue_%s",    par_pdg_str.at(pdg).c_str()), q_nue.at(pdg).c_str());
        t_intrinsic->Draw(Form("nu_e >> h_nuebar_%s", par_pdg_str.at(pdg).c_str()), q_nuebar.at(pdg).c_str());

        h_numu.at(pdg)->SetLineWidth(2);
        h_numu.at(pdg)->SetLineColor(cols.at(pdg));
        h_numubar.at(pdg)->SetLineWidth(2);
        h_numubar.at(pdg)->SetLineColor(cols.at(pdg));
        h_nue.at(pdg)->SetLineWidth(2);
        h_nue.at(pdg)->SetLineColor(cols.at(pdg));
        h_nuebar.at(pdg)->SetLineWidth(2);
        h_nuebar.at(pdg)->SetLineColor(cols.at(pdg));

        if (type == "flugg"){
            h_numu.at(pdg)->Scale(2.6324287);
            h_numubar.at(pdg)->Scale(2.6324287);
            h_nue.at(pdg)->Scale(0.79076896);
            h_nuebar.at(pdg)->Scale(0.79076896);
        }

        // Scale to standard pot
        h_nue.at(pdg)->Scale(0.098239978);
        h_nuebar.at(pdg)->Scale(0.098239978);
        
    }

    leg.at(0)->AddEntry(h_numu.at(7),    Form("%s", par_pdg_str.at(7).c_str()) , "l");
    leg.at(1)->AddEntry(h_numubar.at(7), Form("%s", par_pdg_str.at(7).c_str()) , "l");
    leg.at(2)->AddEntry(h_nue.at(7),     Form("%s", par_pdg_str.at(7).c_str()) , "l");
    leg.at(3)->AddEntry(h_nuebar.at(7),  Form("%s", par_pdg_str.at(7).c_str()) , "l");

    for (unsigned int pdg = 0; pdg < par_pdg.size()-1; pdg++){
        leg.at(0)->AddEntry(h_numu.at(pdg),    Form("%s (%2.1f%%)", par_pdg_str.at(pdg).c_str(), double(100*h_numu.at(pdg)   ->Integral()/h_numu.at(7)   ->Integral() )) , "l");
        leg.at(1)->AddEntry(h_numubar.at(pdg), Form("%s (%2.1f%%)", par_pdg_str.at(pdg).c_str(), double(100*h_numubar.at(pdg)->Integral()/h_numubar.at(7)->Integral() )) , "l");
        leg.at(2)->AddEntry(h_nue.at(pdg),     Form("%s (%2.1f%%)", par_pdg_str.at(pdg).c_str(), double(100*h_nue.at(pdg)    ->Integral()/h_nue.at(7)    ->Integral() )) , "l");
        leg.at(3)->AddEntry(h_nuebar.at(pdg),  Form("%s (%2.1f%%)", par_pdg_str.at(pdg).c_str(), double(100*h_nuebar.at(pdg) ->Integral()/h_nuebar.at(7) ->Integral() )) , "l");
    }

    TCanvas * c;
    gStyle->SetOptStat(0);
    _util.CreateDirectory("Rates/");
    
    for (unsigned int f = 0; f < leg.size(); f++){
        std::cout << "Plotting..." << f << std::endl; 
        
        c = new TCanvas();
        for (unsigned int pdg = 0; pdg < par_pdg.size(); pdg++){
            
            if (f == 0){
                h_numu.at(pdg)->SetMaximum(12000);
                h_numu.at(pdg)->Draw("hist,same");
                _util.IncreaseLabelSize(h_numu.at(pdg), c);
                leg.at(f)->Draw();
            }
            else if (f == 1){
                h_numubar.at(pdg)->SetMaximum(3500);
                h_numubar.at(pdg)->Draw("hist,same");
                _util.IncreaseLabelSize(h_numubar.at(pdg),c);
                leg.at(f)->Draw();
            }
            else if (f==2){
                h_nue.at(pdg)->SetMaximum(900);
                h_nue.at(pdg)->Draw("hist,same");
                _util.IncreaseLabelSize(h_nue.at(pdg),c);
                leg.at(f)->Draw();
            }
            else{
                h_nuebar.at(pdg)->SetMaximum(200);
                h_nuebar.at(pdg)->Draw("hist,same");
                _util.IncreaseLabelSize(h_nuebar.at(pdg), c);
                leg.at(f)->Draw();
            }
        }

        _util.Draw_Nu_Mode(c, 0.14, 0.89, 0.34, 0.96);
        c->Print(Form("plots/run%s/Rates/%s_%s_parent.pdf", _util.run_period, flav_norm.at(f).c_str(), type.c_str() ));
        delete c;
    }

    TFile *f_out = TFile::Open("files/flux_rates.root", "UPDATE");
    f_out->cd();
    for (unsigned int pdg = 0; pdg < par_pdg.size(); pdg++){
        h_numu.at(pdg)   ->Write(Form("numu_%s_%s",   par_pdg_str2.at(pdg).c_str(), type.c_str() ), TObject::kOverwrite);
        h_numubar.at(pdg)->Write(Form("numubar_%s_%s",par_pdg_str2.at(pdg).c_str(), type.c_str() ), TObject::kOverwrite);
        h_nue.at(pdg)    ->Write(Form("nue_%s_%s",    par_pdg_str2.at(pdg).c_str(), type.c_str() ), TObject::kOverwrite);
        h_nuebar.at(pdg) ->Write(Form("nuebar_%s_%s", par_pdg_str2.at(pdg).c_str(), type.c_str() ), TObject::kOverwrite);
    }
    f_out->Close();




}
// -----------------------------------------------------------------------------
void UtilityPlotter::PlotBeamSimRates(){


    TFile *f_in = TFile::Open("files/flux_rates.root", "READ");
    std::vector<std::string> flav_norm = {"numu", "numubar", "nue", "nuebar"};
    std::vector<std::string> par_pdg_str2 = {"pip", "pim", "mup", "mum","kp","km", "kl", "all" };
    std::vector<std::string> flav = {"#nu_{#mu}", "#bar{#nu}_{#mu}", "#nu_{e}", "#bar{#nu}_{e}"};

    // Select which flavour to get from the file
    std::string par = "all";

    std::vector<TH1D*> h_rate_dk2nu(flav_norm.size());
    std::vector<TH1D*> h_rate_flugg(flav_norm.size());
    for (unsigned int f = 0; f < h_rate_dk2nu.size(); f++){
        h_rate_dk2nu.at(f) = (TH1D*)f_in->Get(Form("%s_%s_dk2nu", flav_norm.at(f).c_str(), par.c_str()));
        h_rate_dk2nu.at(f)->SetLineColor(kRed+2);

        h_rate_flugg.at(f) = (TH1D*)f_in->Get(Form("%s_%s_flugg", flav_norm.at(f).c_str(), par.c_str()));
        h_rate_flugg.at(f)->SetLineColor(kGreen+2);

        h_rate_dk2nu.at(f)->SetTitle(Form("%s; Neutrino Energy [GeV]; Entries", flav.at(f).c_str()));
    }

    

    TLegend* leg;
    leg = new TLegend(0.64, 0.65, 0.84, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->AddEntry(h_rate_dk2nu.at(0), "dk2nu", "l");
    leg->AddEntry(h_rate_flugg.at(0), "flugg", "l");

    TCanvas * c;
    gStyle->SetOptStat(0);

    for (unsigned int f = 0; f < h_rate_dk2nu.size(); f++){
        
        c = new TCanvas();
        _util.IncreaseLabelSize(h_rate_dk2nu.at(f),c);
        c->SetLeftMargin(0.17);

        if (f == 0){
            h_rate_dk2nu.at(f)->SetMaximum(12000);
        }
        else if (f == 1){
            h_rate_dk2nu.at(f)->SetMaximum(3500);
        }
        else if (f==2){
            h_rate_dk2nu.at(f)->SetMaximum(900);
        }
        else{
            h_rate_dk2nu.at(f)->SetMaximum(200);
        }

        if (f ==  2|| f == 0) gStyle->SetTitleH(0.1);
	    else gStyle->SetTitleH(0.07);
        
        h_rate_dk2nu.at(f)->Draw("hist");
        h_rate_flugg.at(f)->Draw("hist,same");
        leg->Draw();

        _util.Draw_Nu_Mode(c, 0.14, 0.89, 0.34, 0.96);

        c->Print(Form("plots/run%s/Rates/%s_comparison.pdf", _util.run_period, flav_norm.at(f).c_str() ));
        delete c;
    }

}
// -----------------------------------------------------------------------------
void UtilityPlotter::PrintXSecResults(){

    TFile *f_xsec = TFile::Open("files/xsec_result_run1.root", "READ");


    TH2D* h_temp_2D;
    TH1D* h_temp;

    // Total Covariance Matrix
    h_temp_2D = (TH2D*)f_xsec->Get(Form("%s/wiener/h_data_cov_tot_unfolded",_util.xsec_var));
    TH2D* h_cov_reco = (TH2D*)h_temp_2D->Clone();
    h_cov_reco->SetDirectory(0);

    // Ac
    h_temp_2D = (TH2D*)f_xsec->Get(Form("%s/wiener/h_ac",_util.xsec_var));
    TH2D* h_ac = (TH2D*)h_temp_2D->Clone();
    h_ac->SetDirectory(0);

    // Response
    h_temp_2D = (TH2D*)f_xsec->Get(Form("%s/wiener/h_response",_util.xsec_var));
    TH2D* h_response = (TH2D*)h_temp_2D->Clone();
    h_response->SetDirectory(0);

    // Data XSec
    h_temp  = (TH1D*)f_xsec->Get(Form("%s/wiener/h_data_xsec_unfolded", _util.xsec_var));
    TH1D* unf = (TH1D*)h_temp->Clone();
    unf->SetDirectory(0);
    _util.UndoBinWidthScaling(unf);

    TH1D* unf_width = (TH1D*)h_temp->Clone(); // bin width scaled
    unf_width->SetDirectory(0);

    // Convert the Covariance Matrix-- switching from not binwidth scaled to bin-width scaled
    _util.ConvertCovarianceUnits(h_cov_reco,
                        unf,
                        unf_width);    

    // Load in the histograms for the efficeincy
    TFile *f_eff = TFile::Open(Form("files/nuexsec_mc_run%s.root", _util.run_period));
    TH1D* h_nue_den, *h_nue_num, *h_nuebar_den, *h_nuebar_num;


    if (std::string(_util.xsec_var) == "elec_E"){
        h_nue_den  = (TH1D*)f_eff->Get("TEff/h_true_elec_E_rebin_nue_Unselected");
        h_nue_num  = (TH1D*)f_eff->Get("TEff/h_true_elec_E_rebin_nue_dEdx_max_no_tracks");
        h_nuebar_den  = (TH1D*)f_eff->Get("TEff/h_true_elec_E_rebin_nuebar_Unselected");
        h_nuebar_num  = (TH1D*)f_eff->Get("TEff/h_true_elec_E_rebin_nuebar_dEdx_max_no_tracks");
    }
    else if (std::string(_util.xsec_var) == "elec_cang"){
        h_nue_den  = (TH1D*)f_eff->Get("TEff/h_eff_cosine_beta_rebin_nue_Unselected");
        h_nue_num  = (TH1D*)f_eff->Get("TEff/h_eff_cosine_beta_rebin_nue_dEdx_max_no_tracks");
        h_nuebar_den  = (TH1D*)f_eff->Get("TEff/h_eff_cosine_beta_rebin_nuebar_Unselected");
        h_nuebar_num  = (TH1D*)f_eff->Get("TEff/h_eff_cosine_beta_rebin_nuebar_dEdx_max_no_tracks");
    }

    h_nue_num->Divide(h_nue_den);
    h_nuebar_num->Divide(h_nuebar_den);

    // Print the x-section values
    std::cout << "Printing XSection Values and Uncertainties"<< std::endl;
    std::cout << "Bin & Bin Range & XSec & Err & nue eff & nuebar eff" << std::endl; 
    for (int bin = 1; bin < unf->GetNbinsX()+1; bin++){
        std::cout << bin << std::fixed << std::setprecision(2)<< " & ["  <<  unf_width->GetBinLowEdge(bin) << ", " << unf_width->GetBinLowEdge(bin+1) << ") & " <<  unf_width->GetBinContent(bin) << " & " << std::sqrt(h_cov_reco->GetBinContent(bin, bin)) << " & " << h_nue_num->GetBinContent(bin) << " & " <<h_nuebar_num->GetBinContent(bin) << " \\\\"<< std::endl;
    }

    // Now Print the Covariance Matrix
    std::cout << "\nPrinting Cov matric elements\n\n"<< std::endl;

    for (int j = 1; j < h_cov_reco->GetNbinsX()+1;j++){

        for (int i = 1; i < h_cov_reco->GetNbinsY()+1;i++){
            
            std::cout << std::fixed << std::setprecision(4)<< i << ", " << j <<", " << h_cov_reco->GetBinContent(i,j) << std::endl;
        }

    }

    // Now Print the AC
    std::cout << "\nPrinting AC elements\n\n"<< std::endl;

    for (int j = 1; j < h_ac->GetNbinsX()+1;j++){

        for (int i = 1; i < h_ac->GetNbinsY()+1;i++){
            
            std::cout << std::fixed << std::setprecision(4) << i << ", " << j <<", " << h_ac->GetBinContent(i,j) << std::endl;
        }

    }

    // Now Print the Response
    std::cout << "\nPrinting Response elements\n\n"<< std::endl;

    for (int j = 1; j < h_response->GetNbinsX()+1;j++){

        for (int i = 1; i < h_response->GetNbinsY()+1;i++){
            
            std::cout << std::fixed << std::setprecision(4)<< i << ", " << j <<", " << h_response->GetBinContent(i,j) << std::endl;
        }

    }

    // Write the histograms to file
    TFile *f_result = TFile::Open("files/xsec_result_run1_paper.root", "UPDATE");
    
    if (std::string(_util.xsec_var) == "elec_E"){
        h_cov_reco->Write("unf_cov_energy",TObject::kOverwrite);
        h_ac->Write("ac_energy",TObject::kOverwrite);
        unf_width->Write("unf_data_xsec_energy",TObject::kOverwrite);
        // truexsec->Write("mc_xsec_true_genie_v3_0_6_energy",TObject::kOverwrite);
    }
    else {
        h_cov_reco->Write("unf_cov_angle",TObject::kOverwrite);
        h_ac->Write("ac_angle",TObject::kOverwrite);
        unf_width->Write("unf_data_xsec_angle",TObject::kOverwrite);
        // truexsec->Write("mc_xsec_true_genie_v3_0_6_angle",TObject::kOverwrite);
    }
    
    f_result->Close();

}
// -----------------------------------------------------------------------------
void UtilityPlotter::PrintFluxValues(){

    // Load in the Flux file
    TFile *f_flux = TFile::Open("Systematics/output_fhc_uboone_run0.root","READ");

    TH1D* h_nue    = (TH1D*)f_flux->Get("nue/Detsmear/nue_CV_AV_TPC_5MeV_bin");
    TH1D* h_nuebar = (TH1D*)f_flux->Get("nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin");
    
    h_nue->Rebin(2);
    h_nuebar->Rebin(2);

    // Print the Nue values
    std::cout << "\nPrinting nue flux elements\n\n"<< std::endl;

    for (int bin = 1; bin < h_nue->GetNbinsX()+1; bin++){
        if (h_nue->GetBinLowEdge(bin) >= 10.0)
            continue;

        std::cout << h_nue->GetBinLowEdge(bin) << ", " << h_nue->GetBinLowEdge(bin+1) << ", " << h_nue->GetBinContent(bin) << std::endl;
    }

    // Print the Nuebar values
    std::cout << "\nPrinting nuebar flux elements\n\n"<< std::endl;

    for (int bin = 1; bin < h_nuebar->GetNbinsX()+1; bin++){
        if (h_nuebar->GetBinLowEdge(bin) >= 10.0)
            continue;

        std::cout << h_nuebar->GetBinLowEdge(bin) << ", " << h_nuebar->GetBinLowEdge(bin+1) << ", " << h_nuebar->GetBinContent(bin) << std::endl;
    }

    std::cout <<"\n\n\n" << std::endl;

    // Write the histograms to file
    TFile *f_result = TFile::Open("files/xsec_result_run1_paper.root", "UPDATE");
    h_nue->Write("nue_flux",TObject::kOverwrite);
    h_nuebar->Write("nuebar_flux",TObject::kOverwrite);
    f_result->Close();


    // Make a plot of the nue to nuebar fluxes
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas();
    h_nue->Rebin(10);
    h_nuebar->Rebin(10);
    h_nue->Divide(h_nuebar);
    
    h_nue->SetLineWidth(2);
    h_nue->SetLineColor(kBlue+2);
    h_nue->SetTitle(";Electron Neutrino Energy [GeV]; Ratio #nu_{e}/#bar{#nu}_{e} Flux");
    h_nue->GetXaxis()->SetLabelSize(0.04);
    h_nue->GetXaxis()->SetTitleSize(0.05);
    h_nue->GetYaxis()->SetLabelSize(0.05);
    h_nue->GetYaxis()->SetTitleSize(0.05);
    h_nue->Draw("hist");
    h_nue->GetXaxis()->SetRangeUser(0,5);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);

    _util.Draw_ubooneSim(c, 0.33, 0.925, 0.33, 0.925);

    c->Print("plots/nue_fhc_flux_ratio.pdf");




}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
