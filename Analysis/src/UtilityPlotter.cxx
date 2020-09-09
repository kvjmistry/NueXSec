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
void UtilityPlotter::GetFitResult(double &mean, double &sigma, float bin_lower_edge, float bin_upper_edge, TTree* tree, bool save_hist, bool &converged, bool draw_fit_results){
    
    TCut generic_query = "(classification.c_str()==\"nue_cc\" || classification.c_str()==\"nuebar_cc\") && !gen"; // This gets selected signal events
    TCut bin_query = Form("shr_energy_cali > %f && shr_energy_cali < %f", bin_lower_edge, bin_upper_edge);
    
    TCanvas * c = new TCanvas(Form("c_%f_%f", bin_upper_edge, sigma), "c", 500, 500);

    // Get the histogram from the pad
    TH1D *htemp;
    
    htemp = new TH1D("htemp","", 80, 0, 4.0); // Set the binnning

    // Draw the Query adn put into histogram
    tree->Draw("elec_e >> htemp", generic_query && bin_query);
    
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
        c->Print(Form("plots/run%s/Binning/bins_%0.2fGeV_to_%0.2f_GeV.pdf",_util.run_period, bin_lower_edge, bin_upper_edge ));
    } 

    delete htemp;
    delete c;


}
// -----------------------------------------------------------------------------
void UtilityPlotter::OptimiseBins(){

    // Create the Bins directory for saving the plots to
    _util.CreateDirectory("Binning");

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
                GetFitResult(mean, sigma, 2.32, 4.0, tree, true, converged, false);
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
    if      (variable_str == "reco_e") tree->Draw("(shr_energy_cali - elec_e) / shr_energy_cali >> htemp", generic_query && bin_query);
    else if (variable_str == "true_e") tree->Draw("(shr_energy_cali - elec_e) / elec_e >> htemp", generic_query && bin_query);
    else if (variable_str == "purity") tree->Draw("shr_bkt_purity >> htemp", generic_query && bin_query);
    else if (variable_str == "completeness") tree->Draw("shr_bkt_completeness >> htemp", generic_query && bin_query);
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
    _util.GetFile(f, "Systematics/output_fhc_uboone_run0.root");
    
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

    h_nue->GetXaxis()->SetRangeUser(0,4);
    // h_nue->GetYaxis()->SetRangeUser(0,150.0e6);
    h_nue->SetLineColor(kBlue+2);
    h_nue->SetFillColor(17);
    
    h_nue->GetXaxis()->SetTitleFont(46);
    h_nue->GetXaxis()->SetTitleSize(18);

    if (!draw_averge) h_nue->Draw("hist");

    TH1D* h_nue_clone = (TH1D*)h_nue->Clone("h_nue_clone");
    
    // Nuebar flux
    h_nuebar->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu}_{e} / cm^{2} / 5 MeV / 0.9 #times 10^{20} POT");
    h_nuebar->GetXaxis()->SetRangeUser(0,4);
    h_nuebar->SetLineColor(kGreen+2);
    h_nuebar->SetFillColor(16);
    if (!draw_averge)h_nuebar->Draw("hist,same");

    // Define the threshold 
    double threshold_energy =  0.125; // Current threshold
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
void UtilityPlotter::StudySlipStackingFlashTime(){

    // Load in the TFiles for Run 1 MC, dirt, Data, EXT
    std::vector<TFile*> f_v;
    std::vector<TTree*> tree_v;

    f_v.resize(_util.k_type_MAX);
    tree_v.resize(_util.k_type_MAX);

    // Get the TFiles
    _util.GetFile(f_v.at(_util.k_mc),   "../ntuples/neutrinoselection_filt_run1_overlay_weight.root");
    _util.GetFile(f_v.at(_util.k_dirt), "../ntuples/neutrinoselection_filt_run1_dirt_overlay.root");
    _util.GetFile(f_v.at(_util.k_data), "../ntuples/neutrinoselection_filt_run1_beamon_beamgood.root");
    _util.GetFile(f_v.at(_util.k_ext),  "../ntuples/neutrinoselection_filt_run1_beamoff.root");
    
    // Now get the trees
    for (unsigned int f = 0; f <  f_v.size(); f++){
        f_v.at(f)->cd();
        _util.GetTree(f_v.at(f), tree_v.at(f), "nuselection/NeutrinoSelectionFilter");
    }


    // Redefine the scale factors depending on what POT mode we want to plot
    double Data_POT = 0.9e20;
    // double Data_Trig = 1111111;
    // double mc_scale_factor     = Data_POT  / _util.config_v.at(_util.k_Run1_MC_POT);
    // double ext_scale_factor    = Data_POT  / _util.config_v.at(_util.k_Run1_Dirt_POT);
    // double dirt_scale_factor   = Data_Trig / _util.config_v.at(_util.k_Run1_EXT_trig); 

    

    // Now create the histograms to query
    std::vector<TH1D*> hist;
    hist.resize(_util.k_type_MAX);

    for (unsigned int type = 0; type <  hist.size(); type++){
        hist.at(type) = new TH1D( Form("h_flash_time_%s", _util.type_prefix.at(type).c_str()), "; Flash Time [us]; Entries", 50, 0, 25 );
    }

    TCut q_data = "flash_time > 0"; // choose query to get the relavent flash period
    TCut q_other = "flash_time > 0 && flash_time < 20";

    TCanvas * c = new TCanvas("c", "c", 500, 500);
    
    TH1D* htemp = new TH1D("htemp", "; Flash Time [us]; Entries", 50, 0, 25 );
    f_v.at(_util.k_data)->cd();
    tree_v.at(_util.k_data)->Draw("flash_time >> htemp", q_data);
    htemp = (TH1D*) gDirectory->Get("htemp");
    hist.at(_util.k_data)= (TH1D*)htemp->Clone();

    TH1D* htemp_ext = new TH1D("htemp_ext", "; Flash Time [us]; Entries", 50, 0, 25 );
    f_v.at(_util.k_ext)->cd();
    tree_v.at(_util.k_ext)->Draw("(flash_time - 0.359) >> htemp_ext", q_other);
    htemp_ext = (TH1D*) gDirectory->Get("htemp_ext");
    hist.at(_util.k_ext)= (TH1D*)htemp_ext->Clone();
    // hist.at(_util.k_ext)->Scale(_util.ext_scale_factor);

    TH1D* htemp_mc = new TH1D("htemp_mc", "; Flash Time [us]; Entries", 50, 0, 25 );
    f_v.at(_util.k_mc)->cd();
    tree_v.at(_util.k_mc)->Draw("(flash_time*ppfx_cv*weightSplineTimesTune + 0.055- 0.359) >> htemp_mc", q_other);
    htemp_mc = (TH1D*) gDirectory->Get("htemp_mc");
    hist.at(_util.k_mc)= (TH1D*)htemp_mc->Clone();
    // hist.at(_util.k_mc)->Scale(_util.mc_scale_factor);

     TH1D* htemp_dirt = new TH1D("htemp_dirt", "; Flash Time [us]; Entries", 50, 0, 25 );
    f_v.at(_util.k_dirt)->cd();
    tree_v.at(_util.k_dirt)->Draw("(flash_time*ppfx_cv*weightSplineTimesTune + 0.055 - 0.359) >> htemp_dirt", q_other);
    htemp_dirt = (TH1D*) gDirectory->Get("htemp_dirt");
    hist.at(_util.k_dirt)= (TH1D*)htemp_dirt->Clone();
    // hist.at(_util.k_dirt)->Scale(_util.dirt_scale_factor);


    std::vector<double> hist_integrals(_util.k_type_MAX, 0.0); // The integrals of all the histograms
    double integral_mc_ext = 0.0;

    TH1D *h_ratio;
    TH1D *h_ratio_error;
    TH1D *h_mc_ext_sum;

    TPad *topPad;
    TPad *bottomPad;
    THStack *h_stack = new THStack();

    topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    _util.SetTPadOptions(topPad, bottomPad);

    for (unsigned int i = 0; i < hist.size(); i++)
    {

        if (i == _util.k_data)
        {

            hist.at(i)->SetStats(kFALSE);
            hist.at(_util.k_data)->SetMarkerStyle(20);
            hist.at(_util.k_data)->SetMarkerSize(0.5);
            hist_integrals.at(_util.k_data) = hist.at(_util.k_data)->Integral();
        }

        // Scale EXT
        else if (i == _util.k_ext)
        {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(_util.ext_scale_factor);
            hist.at(_util.k_ext)->SetFillColor(41);
            hist.at(_util.k_ext)->SetFillStyle(3345);
            hist_integrals.at(_util.k_ext) = hist.at(_util.k_ext)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_ext);
        }

        // Scale Dirt
        else if (i == _util.k_dirt)
        {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(_util.dirt_scale_factor);
            hist.at(_util.k_dirt)->SetFillColor(2);
            hist.at(_util.k_dirt)->SetFillStyle(3354);
            hist_integrals.at(_util.k_dirt) = hist.at(_util.k_dirt)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_dirt);
        }

        // Scale MC
        else
        {
            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(_util.mc_scale_factor);
            hist.at(i)->SetFillColor(30);
            hist_integrals.at(_util.k_mc) = hist.at(_util.k_mc)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_mc);
        }
    }

    // Add the histograms to the stack
    h_stack->Add(hist.at(_util.k_ext));
    h_stack->Add(hist.at(_util.k_mc));
    h_stack->Add(hist.at(_util.k_dirt));

    h_stack->Draw("hist");
    hist.at(_util.k_data)->Draw("same PE");
    
    h_stack->GetYaxis()->SetTitle("Entries");

    h_stack->GetYaxis()->SetTitleSize(0.05);
    h_stack->GetYaxis()->SetLabelSize(0.05);
    h_stack->GetXaxis()->SetLabelSize(0);
    h_stack->GetXaxis()->SetRangeUser(0, 23);
    // if (histname == "h_flash_time_single_bin") h_stack->GetXaxis()->SetRangeUser(5.6,15.4);

    // MC error histogram ------------------------------------------------------
    TH1D *h_error_hist = (TH1D *)hist.at(_util.k_mc)->Clone("h_error_hist");

    for (unsigned int i = 0; i < hist.size(); i++)
    {
        if (i == _util.k_data)
            continue; // Dont use the data
        if (i == _util.k_mc)
            continue; // Aleady got this histogram from the clone

        h_error_hist->Add(hist.at(i), 1);
    }

    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2, same");

    TLegend *leg_stack;
    // if (histname == "h_flash_time_single_bin") leg_stack = new TLegend(0.7, 0.2, 0.9, 0.45);
    leg_stack = new TLegend(0.7, 0.6, 0.9, 0.85);
    leg_stack->SetBorderSize(0);
    // if (histname != "h_flash_time_single_bin")leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(hist.at(_util.k_data), "Beam-On Data", "lep");
    leg_stack->AddEntry(hist.at(_util.k_dirt), "Out-of Cryo", "f");
    leg_stack->AddEntry(hist.at(_util.k_mc),   "Overlay", "f");
    leg_stack->AddEntry(hist.at(_util.k_ext),  "Beam-Off Data", "f");

    leg_stack->Draw();

    bottomPad->cd();

    h_ratio = (TH1D *)hist.at(_util.k_data)->Clone("h_ratio");
    h_mc_ext_sum = (TH1D *)hist.at(_util.k_mc)->Clone("h_mc_ext_sum");

    for (unsigned int i = 0; i < hist.size(); i++)
    {
        if (i == _util.k_data || i == _util.k_mc)
            continue; // Dont use the data and nue cc because already been cloned
        h_mc_ext_sum->Add(hist.at(i), 1);
    }

    // h_ratio->Add(h_mc_ext_sum, -1);
    h_ratio->Divide(h_mc_ext_sum);

    h_ratio->GetXaxis()->SetLabelSize(0.13);
    h_ratio->GetXaxis()->SetTitleOffset(0.9);
    h_ratio->GetXaxis()->SetTitleSize(0.13);
    h_ratio->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

    // For percent difference
    // h_ratio->GetYaxis()->SetTitle("(Data - MC) / MC ");
    // h_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);

    // For ratio

    h_ratio->GetYaxis()->SetRangeUser(0.80, 1.20);
    h_ratio->GetXaxis()->SetRangeUser(0, 23);
    // if (histname == "h_flash_time_single_bin") h_ratio->GetXaxis()->SetRangeUser(5.6,15.4);
    h_ratio->GetYaxis()->SetTitle("#frac{Beam-On}{(Overlay + Beam-Off)}");

    h_ratio->GetYaxis()->SetLabelSize(0.13);
    h_ratio->GetYaxis()->SetTitleSize(0.07);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->SetTitle(" ");
    h_ratio->GetYaxis()->CenterTitle();
    h_ratio->Draw("E");

    // Draw the error hist
    h_ratio_error = (TH1D *)h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Draw Data to MC ratio
    _util.Draw_Data_MC_Ratio(c, double(hist_integrals.at(_util.k_data) * 1.0 / integral_mc_ext * 1.0), 0.34, 0.936, 0.34, 0.936);

    // Draw other data specifc quantities
    _util.Draw_Data_POT(c, Data_POT, 0.45, 0.915, 0.45, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);

    // if (area_norm)
    //     Draw_Area_Norm(c);

    // c->Print(print_name);

    c->Print("plots/test.pdf");
    







}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------