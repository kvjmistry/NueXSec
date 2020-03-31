#include "../include/histogram_plotter.h"

// -----------------------------------------------------------------------------
histogram_plotter::~histogram_plotter(){ 
    
    // Make sure the file is closed
    // f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void histogram_plotter::MakeHistograms(const char * hist_file_name, const char *run_period, const std::vector<double> _config, int weight_cfg){

    std::cout << "Creating histograms and making plots" << std::endl;
    
    double Data_POT;

    // Set the scale factors
    if (strcmp(run_period, "1") == 0){
        mc_scale_factor     = _config.at(_util.k_config_Run1_Data_POT)  / _config.at(_util.k_config_Run1_MC_POT) ;
        dirt_scale_factor   = _config.at(_util.k_config_Run1_Data_POT)  / _config.at(_util.k_config_Run1_Dirt_POT);
        intime_scale_factor = _config.at(_util.k_config_Run1_Data_trig) / _config.at(_util.k_config_Run1_EXT_trig);
        Data_POT = _config.at(_util.k_config_Run1_Data_POT); // Define this variable here for easier reading
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
        exit(1);
    }
    

    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << "Scale Factors:\n" <<
    "MC Scale factor:   "   << mc_scale_factor     << "\n" <<
    "Dirt Scale factor: "   << dirt_scale_factor   << "\n" <<
    "EXT Scale factor:  "   << intime_scale_factor << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    Initalise(hist_file_name, run_period, mc_scale_factor, intime_scale_factor, dirt_scale_factor, weight_cfg);

    // Loop over the cuts and plot histograms by plot type
    for (unsigned int i = 0 ; i < _util.k_cuts_MAX; i++){
        
        // Create a set of strings for creating a dynamic directory
        // Directory structure that is created will take the form plots/<cut>/
        std::string a = "if [ ! -d \"plots/";
        std::string b = "run" + std::string(run_period) + "/cuts/" + _util.cut_dirs.at(i);
        std::string c = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
        std::string d = "run" + std::string(run_period) + "/cuts/" + _util.cut_dirs.at(i);
        std::string e = "; fi";
        std::string command = a + b + c + d + e ;
        system(command.c_str()); 

        // Call the Make stack function for all the plots we want
        // CallMakeStack(run_period, i, Data_POT);
        
    }

    // Create a set of strings for creating a dynamic directory
    // Directory structure that is created will take the form plots/<cut>/
    std::string a = "if [ ! -d \"plots/";
    std::string b = "run" + std::string(run_period) + "/" + "Flash";
    std::string c = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
    std::string d = "run" + std::string(run_period) + "/" + "Flash";
    std::string e = "; fi";
    std::string command = a + b + c + d + e ;
    system(command.c_str()); 

    // Flash time plots
    MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_time.pdf", run_period), "h_flash_time");
    MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_time_sid1.pdf", run_period), "h_flash_time_sid1"); // Slice ID 1
    MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_time_sid0.pdf", run_period), "h_flash_time_sid0"); // Slice ID 0

    // Flash PE plots
    MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_pe.pdf", run_period), "h_flash_pe");
    MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_pe_sid1.pdf", run_period), "h_flash_pe_sid1"); // Slice ID 1
    MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_pe_sid0.pdf", run_period), "h_flash_pe_sid0"); // Slice ID 0

    // On beam minus off beam plots

    // Flash time
    MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_time_OMO.pdf", run_period), "h_flash_time");
    MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_time_OMO_sid1.pdf", run_period), "h_flash_time_sid1"); // Slice ID 1
    MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_time_OMO_sid0.pdf", run_period), "h_flash_time_sid0"); // Slice ID 0

    // Flash PE plots
    MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_pe_OMO.pdf", run_period), "h_flash_pe");
    MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_pe_OMO_sid1.pdf", run_period), "h_flash_pe_sid1"); // Slice ID 1
    MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_pe_OMO_sid0.pdf", run_period), "h_flash_pe_sid0"); // Slice ID 0

    a = "if [ ! -d \"plots/";
    b = "run" + std::string(run_period) + "/" + "Efficiency";
    c = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
    d = "run" + std::string(run_period) + "/" + "Efficiency";
    e = "; fi";
    command = a + b + c + d + e ;
    system(command.c_str()); 


    MakeEfficiencyPlot(Form("plots/run%s/Efficiency/Integrated_Efficiency_Purity.pdf", run_period));

    MakeEfficiencyPlotByCut(Form("plots/run%s/Efficiency/TEff.pdf", run_period), run_period);


}
// -----------------------------------------------------------------------------
void histogram_plotter::Initalise(const char * hist_file_name, const char *_run_period, double _mc_scale_factor, double _intime_scale_factor, double _dirt_scale_factor, int weight_cfg){ 
    
    std::cout << "Initalising Histogram Plotter..." << std::endl;

    // Set the run period
    run_period = _run_period;

    mc_scale_factor = _mc_scale_factor;
    intime_scale_factor = _intime_scale_factor;
    dirt_scale_factor = _dirt_scale_factor;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject(hist_file_name) ) {
        f_nuexsec = TFile::Open(hist_file_name);
    }
    else {
        std::cout << "Can't find histogram file!! "<<  __PRETTY_FUNCTION__ << std::endl;
        exit(1);
    }

    // Set the weight settings
    if (weight_cfg == 0){
        weight_tune = false;
        weight_ppfx = false;
    }
    else if (weight_cfg == 1){
        weight_tune = true;
        weight_ppfx = true;
    }
    else if (weight_cfg == 2){
        weight_tune = true;
        weight_ppfx = false;
    }
    else if (weight_cfg == 3){
        weight_tune = false;
        weight_ppfx = true;
    }
    else {
        std::cout << "Unknown weight setting specified, using defaults" << std::endl;
    }
}
// -----------------------------------------------------------------------------
std::vector<double> histogram_plotter::Chi2Calc(TH1D * h_mc_ext, TH1D * h_data, const bool area_norm, const double return_norm){
    const int n_bins = h_mc_ext->GetNbinsX();

    const double f_1 = h_mc_ext->Integral();
    const double f_2 = h_data  ->Integral();

    //area normalised?
    TH1D * h_mc_ext_clone = (TH1D*)h_mc_ext->Clone("h_mc_ext_clone");
    TH1D * h_data_clone   = (TH1D*)h_data  ->Clone("h_data_clone");
    
    if(!area_norm) {h_mc_ext_clone->Scale(f_2/f_1); }
    
    if(area_norm) {
        //this keeps them area normalised,
        //but at the original values, not 0->1
        //which messes with the chi2 and p calc
        h_mc_ext_clone->Scale(return_norm);
        h_data_clone  ->Scale(return_norm);

        const double f_1_adj = h_mc_ext->Integral();
        const double f_2_adj = h_data->Integral();
        h_mc_ext_clone->Scale(f_2_adj/f_1_adj);

    }
    
    //h_data_clone->Scale(1./f_2);

    std::vector <double> chi2;
    double chi2_val     = 0;
    double n_mc_ext_val = 0;
    double n_data_val   = 0;
    
    // Loop over each bin
    for ( int i = 1; i < n_bins; i++) {
        const double n_mc_ext = h_mc_ext_clone->GetBinContent(i);
        const double n_data   = h_data_clone  ->GetBinContent(i);

        //don't calculate chi2 for bins where no comparison possible
        if(n_data == 0 || n_mc_ext == 0) { continue; }

        //chi2_val += (pow((n_mc_ext - n_data),2) / n_mc_ext);
        //chi2_val += (pow((n_data - n_mc_ext),2)) / (((n_data * f_2) / pow(f_2, 2)) + ((n_mc_ext * f_1) / pow(f_1, 2)));
        chi2_val += 2 * (n_mc_ext - n_data + (n_data * TMath::Log(n_data/n_mc_ext)));

        n_mc_ext_val += n_mc_ext;
        n_data_val   += n_data;
    }
    
    const double reduced_chi2 = chi2_val / (n_bins - 1);
    const double p = TMath::Prob(chi2_val, n_bins);

    chi2.push_back(reduced_chi2);
    //correct this value back to the un-normalised
    //chi2.push_back(n_mc_ext_val * f_1);
    //chi2.push_back(n_data_val * f_2);
    chi2.push_back(n_mc_ext_val * (f_1 / f_2));
    chi2.push_back(n_data_val);
    chi2.push_back(n_bins - 1);
    chi2.push_back(p);
    return chi2;
}
// -----------------------------------------------------------------------------
void histogram_plotter::Draw_Run_Period(TCanvas* c){
    c->cd();

    TPaveText *pt;

    if (run_period == "1"){
        pt = new TPaveText(0.66, 0.89, 0.86, 0.96,"NDC");
        pt->AddText("Run1");
        pt->SetTextColor(kRed+2);
    }
    else if (run_period == "3b"){
        pt = new TPaveText(0.66, 0.89, 0.86, 0.96,"NDC");
        pt->AddText("Run3b");
        pt->SetTextColor(kBlue+2);
    }
    else {
        pt = new TPaveText(0.66, 0.89, 0.86, 0.96,"NDC");
        pt->AddText("RunXXX");
        pt->SetTextColor(kGreen+2);
    }
    
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.04);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void histogram_plotter::Draw_Data_MC_Ratio(TCanvas* c, double ratio){
    c->cd();

    TPaveText *pt;

    pt = new TPaveText(0.16, 0.88, 0.3, 0.95,"NDC");
    pt->AddText(Form("Data/MC Ratio: %2.2f", ratio));
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void histogram_plotter::Draw_Data_POT(TCanvas* c, double pot){
    c->cd();

    TPaveText *pt;

    // Change scale of POT
    double POT = pot/1.0e20;

    pt = new TPaveText(0.849, 0.485, 0.919, 0.485,"NDC");
    pt->AddText(Form("Data POT: %2.1fe20", POT));
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void histogram_plotter::Draw_WeightLabels(TCanvas* c){
    c->cd();

    TPaveText *pt, *pt2;

    pt = new TPaveText(0.840, 0.44, 0.91, 0.44,"NDC");
    pt->AddText(Form("#muB Genie Tune"));
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    if (weight_tune) pt->Draw();

    pt2 = new TPaveText(0.839, 0.40, 0.909, 0.40,"NDC");
    pt2->AddText(Form("PPFX CV Corr."));
    pt2->SetBorderSize(0);
    pt2->SetFillColor(0);
    pt2->SetFillStyle(0);
    pt2->SetTextSize(0.03);
    if (weight_ppfx) pt2->Draw();

}
// -----------------------------------------------------------------------------
void histogram_plotter::Draw_Area_Norm(TCanvas* c){
    c->cd();

    TPaveText *pt;

    pt = new TPaveText(0.4, 0.916, 0.4, 0.916,"NDC");
    pt->AddText("Area Normalised");
    pt->SetTextColor(kGreen+2);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void histogram_plotter::SetTPadOptions(TPad * topPad, TPad * bottomPad ){

    topPad   ->SetBottomMargin(0.05);
    topPad   ->SetTopMargin(0.15);
    bottomPad->SetTopMargin(0.04);
    bottomPad->SetBottomMargin(0.25);
    bottomPad->SetGridy();
    topPad->SetLeftMargin(0.15);
    topPad->SetRightMargin(0.20 );
    bottomPad->SetLeftMargin(0.15);
    bottomPad->SetRightMargin(0.20 );
    topPad   ->Draw();
    bottomPad->Draw();
    topPad   ->cd();

}
// -----------------------------------------------------------------------------
void histogram_plotter::MakeStack(std::string hist_name, std::string cut_name, bool area_norm, bool logy, double y_scale_factor, const char* x_axis_name,
                                     const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, double Data_POT, const char* print_name ){

    std::vector<TH1D*>  hist(_util.k_classifications_MAX);                      // The vector of histograms from the file for the plot
    std::vector<double> hist_integrals(_util.k_classifications_MAX,0.0);        // The integrals of all the histograms
    double integral_mc_ext{0.0};                                                // Integral of MC + EXT -- needs to be removed
    double y_maximum{0};                                                        // y max for scaling histogram scale
    
    std::vector <double> chi2;
    TPaveText * pt_bottom;

    TH1D * h_ratio;
    TH1D * h_mc_ext_sum;


    // bools for checking if plots exist in the file
    bool found_data = true;
    bool found_ext  = true;
    bool found_dirt = true;
    

    const bool p_value = false; // Choose whether to use a pvalue

    TPad * topPad;
    TPad * bottomPad;
    TCanvas * c       = new TCanvas();
    THStack * h_stack = new THStack();

    // Loop over the classifications and get the histograms
    for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

        // Data
        if (i == _util.k_leg_data){
            
            _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            if (hist.at(i) == NULL){
                found_data = false;
            }
        } 
        // EXT
        else if (i == _util.k_leg_ext){
            
            _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s",  cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            if (hist.at(i) == NULL){
                found_ext = false;
            } 
        }
        // Dirt
        else if (i == _util.k_leg_dirt){
            
            _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            
            if (hist.at(i) == NULL){
                found_dirt = false;
            } 
        }
        // MC
        else {
        
            // MC
            if (hist.at(i) != NULL && ( i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt)) continue;

            _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));

            // Must have MC for this to work for now...
            if (hist.at(i) == NULL) return;
        }
    
    }
    
    // If there is data and ext, then we make the ratio plot too
    if (found_data && found_ext){
        topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        
        SetTPadOptions(topPad, bottomPad );
        
    }
    // Otherwise just use an unsplit canvas
    else {
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.20 );
    }
    
    // Scaling and setting stats
    for (unsigned int i=0; i < hist.size(); i++){
        
        if (i == _util.k_leg_data){
            if (found_data){
                hist.at(i)->SetStats(kFALSE);
                hist_integrals.at(i) = hist.at(i)->Integral();
            } 
        } 
        
        // Scale EXT
        else if (i == _util.k_leg_ext){
            if (found_ext) {

                hist.at(i)->SetStats(kFALSE);
                hist.at(i)->Scale(intime_scale_factor);
                hist_integrals.at(i) = hist.at(i)->Integral();
                integral_mc_ext      += hist.at(i)->Integral();
            }
        }

        // Scale Dirt
        else if (i == _util.k_leg_dirt){
            if (found_dirt) {

                hist.at(i)->SetStats(kFALSE);
                hist.at(i)->Scale(dirt_scale_factor);
                hist_integrals.at(i) = hist.at(i)->Integral();
                integral_mc_ext      += hist.at(i)->Integral();
            }
        }
        
        // Scale MC
        else {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(mc_scale_factor);
            hist_integrals.at(i) = hist.at(i)->Integral();
            integral_mc_ext      += hist.at(i)->Integral();
            // std::cout <<  hist.at(i)->Integral()<< std::endl;
        }

    }

    // Customise the histograms -- put in a function
    hist.at(_util.k_nue_cc)       ->SetFillColor(30);
    hist.at(_util.k_nue_cc_mixed) ->SetFillColor(38);
    hist.at(_util.k_numu_cc)      ->SetFillColor(28);
    hist.at(_util.k_numu_cc_pi0)  ->SetFillColor(42);
    hist.at(_util.k_nc_pi0)       ->SetFillColor(36);
    hist.at(_util.k_cosmic)       ->SetFillColor(1);
    hist.at(_util.k_nc)           ->SetFillColor(46);
    hist.at(_util.k_nu_out_fv)    ->SetFillColor(kViolet-7);
    hist.at(_util.k_numu_cc_pi0)  ->SetFillColor(42);
    hist.at(_util.k_unmatched)    ->SetFillColor(12);
    
    if (found_data){
        hist.at(_util.k_leg_data)     ->SetMarkerStyle(20);
        hist.at(_util.k_leg_data)     ->SetMarkerSize(0.5);
    }
    if (found_ext){
        hist.at(_util.k_leg_ext)      ->SetFillColor(41);
        hist.at(_util.k_leg_ext)      ->SetFillStyle(3345);
    }
    
    if (found_dirt){
        hist.at(_util.k_leg_dirt)     ->SetFillColor(2);
        hist.at(_util.k_leg_dirt)     ->SetFillStyle(3354);
    }
    
    // Normalisation by area
    if (area_norm && found_data) {
        
        if (integral_mc_ext != 0) {
            
            for (unsigned int i=0; i < hist.size(); i++){
                if (i == _util.k_leg_data) continue; // Dont scale the data 
                hist.at(i)->Scale( hist_integrals.at(_util.k_leg_data) / integral_mc_ext );
            }

        }

    }

    // Add the histograms to the stack
    for (unsigned int i=0; i < hist.size(); i++){
        if (i == _util.k_leg_data) continue; // Dont use the data
        if (i == _util.k_leg_ext && !found_ext) continue;   // Skip ext if not there
        if (i == _util.k_leg_dirt && !found_dirt) continue; // skip dirt if not there

        h_stack->Add(hist.at(i));
    }


    // Get the maximum y value for scaling
    if (found_data) y_maximum = std::max(hist.at(_util.k_leg_data)->GetMaximum(), h_stack->GetMaximum());
    else y_maximum = h_stack->GetMaximum();

    // Set the axis in the case of a log plot
    if (logy == true && found_data){
        
        if (hist.at(_util.k_nue_cc)->GetMinimum() != 0.0) hist.at(_util.k_leg_data)->SetMinimum(hist.at(_util.k_nue_cc)->GetMinimum() / 2.); 
        
        if (hist.at(_util.k_nue_cc)->GetMinimum() == 0.0) hist.at(_util.k_leg_data)->SetMinimum(hist.at(_util.k_nue_cc)->GetMinimum() + 0.0001 / 2.); 
        
        hist.at(_util.k_leg_data)->SetMaximum(y_maximum * (y_scale_factor * 500));
        

        h_stack->Draw("hist");

    }
    else if(logy && !found_data) {
        h_stack->SetMinimum(0.1);
        c->SetLogy();
        h_stack->Draw("hist");
    }
    else { // Set the axis in the case of a non-log plot
        h_stack->SetMinimum(0);
        h_stack->SetMaximum(y_maximum * y_scale_factor);
        h_stack->Draw("hist");
    }

    // Set the y axis of the stack
    if(!area_norm) h_stack->GetYaxis()->SetTitle("Entries");
    else           h_stack->GetYaxis()->SetTitle("Entries [A.U.]");
    
    // Customise the stacked histogram
    h_stack->GetYaxis()->SetTitleFont(45);
    h_stack->GetYaxis()->SetTitleSize(18);
    h_stack->GetYaxis()->SetTitleOffset(1.30);
    
    if (found_data && found_ext) h_stack->GetXaxis()->SetLabelOffset(10);
    else h_stack->GetXaxis()->SetTitle(x_axis_name);
    
    if (found_data) hist.at(_util.k_leg_data)->Draw("same PE");

    // MC error histogram ------------------------------------------------------
    TH1D * h_error_hist = (TH1D*) hist.at(_util.k_nue_cc)->Clone("h_error_hist");

    for (unsigned int i=0; i < hist.size(); i++){
        if (i == _util.k_leg_data || i == _util.k_nue_cc) continue; // Dont use the data
        if (i == _util.k_leg_ext && !found_ext) continue;           // Skip ext if not there
        if (i == _util.k_leg_dirt && !found_dirt) continue;         // skip dirt if not there
        
        if (i == _util.k_nue_cc) continue; // Aleady got this histogram from the clone
        
        h_error_hist->Add(hist.at(i), 1);
    }
    
    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2 hist same");

    // Set the legend ----------------------------------------------------------
    TLegend *leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    if (found_data) leg_stack->AddEntry(hist.at(_util.k_leg_data), Form("Data (%2.1f)",                hist_integrals.at(_util.k_leg_data)),     "lep");
    if (found_dirt) leg_stack->AddEntry(hist.at(_util.k_leg_dirt), Form("Dirt (%2.1f)",                hist_integrals.at(_util.k_leg_dirt)),     "f");
    if (found_ext)  leg_stack->AddEntry(hist.at(_util.k_leg_ext),  Form("InTime (EXT) (%2.1f)",        hist_integrals.at(_util.k_leg_ext)),      "f");
    leg_stack->AddEntry(hist.at(_util.k_unmatched),                Form("Unmatched (%2.1f)",           hist_integrals.at(_util.k_unmatched)),    "f");
    leg_stack->AddEntry(hist.at(_util.k_nc_pi0),                   Form("NC #pi^{0} (%2.1f)",          hist_integrals.at(_util.k_nc_pi0)),       "f");
    leg_stack->AddEntry(hist.at(_util.k_nc),                       Form("NC (%2.1f)",                  hist_integrals.at(_util.k_nc)),           "f");
    leg_stack->AddEntry(hist.at(_util.k_numu_cc_pi0),              Form("#nu_{#mu} CC #pi^{0} (%2.1f)",hist_integrals.at(_util.k_numu_cc_pi0)),  "f");
    leg_stack->AddEntry(hist.at(_util.k_numu_cc),                  Form("#nu_{#mu} CC (%2.1f)",        hist_integrals.at(_util.k_numu_cc)),      "f");
    leg_stack->AddEntry(hist.at(_util.k_cosmic),                   Form("Cosmic (%2.1f)",              hist_integrals.at(_util.k_cosmic)),       "f");
    leg_stack->AddEntry(hist.at(_util.k_nu_out_fv),                Form("#nu OutFV (%2.1f)",           hist_integrals.at(_util.k_nu_out_fv)),    "f");
    // leg_stack->AddEntry(hist.at(_util.k_nue_cc_mixed),             Form("#nu_{e} CC Mixed (%2.1f)",    hist_integrals.at(_util.k_nue_cc_mixed)), "f"); // This isnt filled anymore
    leg_stack->AddEntry(hist.at(_util.k_nue_cc),                   Form("#nu_{e} CC (%2.1f)",          hist_integrals.at(_util.k_nue_cc)),       "f");

    leg_stack->Draw();

    if (!logy) h_stack->GetYaxis()->SetRangeUser(0, y_maximum * y_scale_factor);
    else if (logy && found_data)    topPad->SetLogy();
    
    // Calculate the chi2 
    if (found_data) {
        TH1D * h_last = (TH1D*) h_stack->GetStack()->Last();
        chi2  = Chi2Calc(h_last, hist.at(_util.k_leg_data), area_norm, hist_integrals.at(_util.k_data));
    }
    
    // Now create the ratio of data to MC ----------------------------------
    if (found_data) {
        
        bottomPad->cd();

        h_ratio    = (TH1D*) hist.at(_util.k_leg_data)->Clone("h_ratio");
        h_mc_ext_sum = (TH1D*) hist.at(_util.k_nue_cc)  ->Clone("h_mc_ext_sum");
        
        for (unsigned int i=0; i < hist.size(); i++){
            if (i == _util.k_leg_data || i == _util.k_nue_cc ) continue; // Dont use the data and nue cc because already been cloned
            h_mc_ext_sum->Add(hist.at(i), 1);
        }

        h_ratio->Add(h_mc_ext_sum, -1);
        h_ratio->Divide(h_mc_ext_sum);
    
        h_ratio->GetXaxis()->SetLabelSize(12);
        h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_ratio->GetYaxis()->SetLabelSize(11);
        h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
        h_ratio->GetXaxis()->SetTitleOffset(3.0);
        h_ratio->GetXaxis()->SetTitleSize(17);
        h_ratio->GetXaxis()->SetTitleFont(46);
        h_ratio->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

        h_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);
        h_ratio->GetXaxis()->SetTitle(x_axis_name);
        h_ratio->GetYaxis()->SetTitle("(Data - MC) / MC ");
        h_ratio->GetYaxis()->SetTitleSize(13);
        h_ratio->GetYaxis()->SetTitleFont(44);
        h_ratio->GetYaxis()->SetTitleOffset(1.5);
        h_ratio->SetTitle(" ");
        h_ratio->Draw("E");

        // Now doing this stuff on the bottom pad
        //x_min, y_min, x_max, y_max
        // Reduced chi2
        pt_bottom = new TPaveText(.12, .80, .30, .96, "NBNDC");
        std::ostringstream o_string_bottom;
        o_string_bottom.precision(3);
        o_string_bottom << std::fixed;
        o_string_bottom << float(chi2.at(0) * chi2.at(3));
        std::string convert_string_bottom = o_string_bottom.str();

        std::ostringstream o_string3_bottom;
        o_string3_bottom << int(chi2.at(3));
        std::string convert_string3_bottom = o_string3_bottom.str();

        std::string chi2_string_bottom = "#chi_{Stat}^{2}/DOF=(" + convert_string_bottom + "/" + convert_string3_bottom + ")";
        pt_bottom->AddText(chi2_string_bottom.c_str());
        pt_bottom->SetFillStyle(0);
        pt_bottom->SetBorderSize(0);
        // pt_bottom->Draw();
    }

    // Draw the run period on the plot
    Draw_Run_Period(c);
    
    // Draw other data specifc quantities
    if (found_data){
        Draw_Data_MC_Ratio(c, hist_integrals.at(_util.k_leg_data)/integral_mc_ext );
        Draw_Data_POT(c, Data_POT);
        if (area_norm) Draw_Area_Norm(c);
        // Add the weight labels
        Draw_WeightLabels(c);
    }
    
    c->Print(print_name);


}
// -----------------------------------------------------------------------------
void histogram_plotter::CallMakeStack(const char *run_period, int cut_index, double Data_POT){

    // MakeStack(std::string hist_name, std::string cut_name, bool area_norm, bool logy, const char* x_axis_name, double y_scale_factor, 
    //                             const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, const char* print_name )

    bool area_norm = false;

    // Reco X
    MakeStack("h_reco_vtx_x",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Reco Vertex X [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_vtx_x.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );
    
    // Reco Y
    MakeStack("h_reco_vtx_y",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Reco Vertex Y [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_vtx_y.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Reco Z
    MakeStack("h_reco_vtx_z",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Reco Vertex Z [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_vtx_z.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Reco X SCE
    MakeStack("h_reco_vtx_x_sce",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Reco Vertex X (Space Charge Corr) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_vtx_x_sce.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );
    
    // Reco Y SCE
    MakeStack("h_reco_vtx_y_sce",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Reco Vertex Y (Space Charge Corr) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_vtx_y_sce.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Reco Z SCE
    MakeStack("h_reco_vtx_z_sce",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Reco Vertex Z (Space Charge Corr) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_vtx_z_sce.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );


    // dEdx
    MakeStack("h_reco_dEdx_y_plane",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Collection Plane dEdx (uncalibrated) [MeV/cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_dEdx_y_plane.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // dEdx cali
    MakeStack("h_reco_dEdx_cali_y_plane",_util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Collection Plane dEdx (calibrated) [MeV/cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_dEdx_cali_y_plane.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading Shower Momentum
    MakeStack("h_reco_leading_mom", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm, false, 1.0, "Leading Shower Momentum [MeV/c]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_mom.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // 2D distance largest flash to reco nu vertex
    MakeStack("h_reco_flash_to_vtx_dist", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Flash to Vertex Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_flash_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // 2D distance shower vertex to reco nu vertex
    MakeStack("h_reco_shower_to_vtx_dist", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Shower to Vertex Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_shower_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // 2D distance track vertex to reco nu vertex
    MakeStack("h_reco_trac_util.k_to_vtx_dist", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Track to Vertex Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_trac_util.k_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading Shower hits in all planes
    MakeStack("h_reco_leading_shower_hits_all_planes", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Leading Shower Hits All Planes", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_shower_hits_all_planes.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading Shower hits in collection
    MakeStack("h_reco_leading_shower_hits_collection_plane", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Leading Shower Hits Collection Plane", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_shower_hits_collection_plane.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading Shower opening angle
    MakeStack("h_reco_leading_shower_open_angle", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Leading Shower Open Angle [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_shower_open_angle.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Secondary shower to vertex distance (for events with more than 1 shower)
    MakeStack("h_reco_secondary_shower_to_vtx_dist", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Secondary Shower to Vertex Distance (>1 shower) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_secondary_shower_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading Shower hits per length
    MakeStack("h_reco_leading_shower_hits_per_length", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Leading Shower Hits / Length [cm^{-1}]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_shower_hits_per_length.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Longest track to leading shower length
    MakeStack("h_reco_longest_trac_util.k_leading_shower_length", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Longest Track Length / Leading Shower Length", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_longest_track_leading_shower_length.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Track Containment
    MakeStack("h_reco_track_contained", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Contained Tracks", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_track_contained.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading shower phi
    MakeStack("h_reco_leading_shower_phi", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Leading Shower Phi [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_shower_phi.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading shower theta
    MakeStack("h_reco_leading_shower_theta", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Leading Shower Theta [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_shower_theta.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading shower cos theta
    MakeStack("h_reco_leading_shower_cos_theta", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Leading Shower Cos(#theta)", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_leading_shower_cos_theta.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading shower multiplicity
    MakeStack("h_reco_shower_multiplicity", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Shower Multiplicty", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_shower_multiplicity.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Leading track multiplicity
    MakeStack("h_reco_track_multiplicity", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Track Multiplicty", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_track_multiplicity.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Topological Score
    MakeStack("h_reco_topological_score", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  true, 1.0, "Topological Score", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_topological_score.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Track shower dist
    MakeStack("h_reco_track_shower_dist", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Longest Track Leading Shower Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_track_shower_dist.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );
    
    // Track shower angle
    MakeStack("h_reco_track_shower_angle", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  true, 1.0, "Longest Track Leading Shower Angle [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_track_shower_angle.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Ratio hits from showers to slice
    MakeStack("h_reco_hits_ratio", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Hit Ratio of all Showers and the Slice", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_hits_ratio.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );
    
        // Shower score
    MakeStack("h_reco_shower_score", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Shower Score", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_shower_score.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Track score
    MakeStack("h_reco_track_score", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Track Score", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_track_score.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Calibrated energy of all the showers
    MakeStack("h_reco_shower_energy_tot_cali", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Total Calibrated Energy of all Showers [GeV]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_shower_energy_tot_cali.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    // Total number of hits for the leading showe
    MakeStack("h_reco_shower_hits", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Total Num of hits for the leading Shower", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_shower_hits.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );

    
    // Total number of hits for the leading shower in the collection plane
    MakeStack("h_reco_shower_hits_y_plane", _util.cut_dirs.at(cut_index).c_str(),
                        area_norm,  false, 1.0, "Total Num of hits for the leading Shower in Collection Plane", 0.8, 0.98, 0.87, 0.32, Data_POT,
                        Form("plots/run%s/cuts/%s/reco_shower_hits_y_plane.pdf", run_period, _util.cut_dirs.at(cut_index).c_str()) );



}
// -----------------------------------------------------------------------------
void histogram_plotter::MakeFlashPlot(double Data_POT, const char* print_name, std::string histname){

    std::vector<TH1D*>  hist(_util.k_type_MAX);
    std::vector<double> hist_integrals(_util.k_type_MAX,0.0);        // The integrals of all the histograms
    double integral_mc_ext = 0.0;
    
    TH1D * h_ratio;
    TH1D * h_mc_ext_sum;

    TPad * topPad;
    TPad * bottomPad;
    TCanvas * c       = new TCanvas();
    THStack * h_stack = new THStack();

    for (unsigned int k = 0; k < _util.type_prefix.size(); k++){
        _util.GetHist(f_nuexsec, hist.at(k), Form("Flash/%s_%s", histname.c_str(),_util.type_prefix.at(k).c_str() ));
        if (hist.at(k) == NULL){
            std::cout << "Couldn't get all the flash histograms so exiting function..."<< std::endl;
            return;
        }
    }

    topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    
    SetTPadOptions(topPad, bottomPad );
        
    for (unsigned int i=0; i < hist.size(); i++){
        
        if (i == _util.k_data){
            
            hist.at(i)->SetStats(kFALSE);
            hist.at(_util.k_data)     ->SetMarkerStyle(20);
            hist.at(_util.k_data)     ->SetMarkerSize(0.5);
            hist_integrals.at(_util.k_data) = hist.at(_util.k_data)->Integral();
        } 
        
        // Scale EXT
        else if (i == _util.k_ext){
            
            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(intime_scale_factor);
            hist.at(_util.k_ext)      ->SetFillColor(41);
            hist.at(_util.k_ext)      ->SetFillStyle(3345);
            hist_integrals.at(_util.k_ext) = hist.at(_util.k_ext)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_ext);

        }

        // Scale Dirt
        else if (i == _util.k_dirt){

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(dirt_scale_factor);
            hist.at(_util.k_dirt)     ->SetFillColor(2);
            hist.at(_util.k_dirt)     ->SetFillStyle(3354);
            hist_integrals.at(_util.k_dirt) = hist.at(_util.k_dirt)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_dirt);

        }
        
        // Scale MC
        else {
            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(mc_scale_factor);
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

    // MC error histogram ------------------------------------------------------
    TH1D * h_error_hist = (TH1D*) hist.at(_util.k_mc)->Clone("h_error_hist");

    for (unsigned int i=0; i < hist.size(); i++){
        if (i == _util.k_data) continue; // Dont use the data
        if (i == _util.k_mc) continue;   // Aleady got this histogram from the clone
        
        h_error_hist->Add(hist.at(i), 1);
    }
    
    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2 hist same");    

    TLegend *leg_stack = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(hist.at(_util.k_data), "Data",   "lep");
    leg_stack->AddEntry(hist.at(_util.k_dirt), "Dirt",   "f");
    leg_stack->AddEntry(hist.at(_util.k_mc),   "Overlay","f");
    leg_stack->AddEntry(hist.at(_util.k_ext),  "InTime", "f");


    leg_stack->Draw();

    bottomPad->cd();

    

    h_ratio      = (TH1D*) hist.at(_util.k_data)->Clone("h_ratio");
    h_mc_ext_sum = (TH1D*) hist.at(_util.k_mc)  ->Clone("h_mc_ext_sum");
    
    for (unsigned int i=0; i < hist.size(); i++){
        if (i == _util.k_data || i == _util.k_mc ) continue; // Dont use the data and nue cc because already been cloned
        h_mc_ext_sum->Add(hist.at(i), 1);
    }

    h_ratio->Add(h_mc_ext_sum, -1);
    h_ratio->Divide(h_mc_ext_sum);

    h_ratio->GetXaxis()->SetLabelSize(12);
    h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetYaxis()->SetLabelSize(11);
    h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetXaxis()->SetTitleOffset(3.0);
    h_ratio->GetXaxis()->SetTitleSize(17);
    h_ratio->GetXaxis()->SetTitleFont(46);
    h_ratio->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

    h_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);
    h_ratio->GetYaxis()->SetTitle("(Data - MC) / MC ");
    h_ratio->GetYaxis()->SetTitleSize(13);
    h_ratio->GetYaxis()->SetTitleFont(44);
    h_ratio->GetYaxis()->SetTitleOffset(1.5);
    h_ratio->SetTitle(" ");
    h_ratio->Draw("E");

    // Draw the run period on the plot
    Draw_Run_Period(c);

    // Draw Data to MC ratio
    Draw_Data_MC_Ratio(c, double(hist_integrals.at(_util.k_data)*1.0/integral_mc_ext*1.0) );
    
    // Draw other data specifc quantities
    Draw_Data_POT(c, Data_POT);

    // Add the weight labels
    Draw_WeightLabels(c);
    
    c->Print(print_name);

}
// -----------------------------------------------------------------------------
void histogram_plotter::MakeFlashPlotOMO(double Data_POT, const char* print_name, std::string histname){

    std::vector<TH1D*>  hist(_util.k_type_MAX);
    std::vector<double> hist_integrals(_util.k_type_MAX,0.0);        // The integrals of all the histograms
    double integral_mc_ext = 0.0;
    
    TH1D * h_ratio;
    TH1D * h_mc_ext_sum;

    TPad * topPad;
    TPad * bottomPad;
    TCanvas * c       = new TCanvas();
    THStack * h_stack = new THStack();

    for (unsigned int k = 0; k < _util.type_prefix.size(); k++){
        _util.GetHist(f_nuexsec, hist.at(k), Form("Flash/%s_%s", histname.c_str(),_util.type_prefix.at(k).c_str() ));
        if (hist.at(k) == NULL){
            std::cout << "Couldn't get all the flash histograms so exiting function..."<< std::endl;
            return;
        }
    }

    topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    
    SetTPadOptions(topPad, bottomPad );
        
    for (unsigned int i=0; i < hist.size(); i++){
        
        if (i == _util.k_data){
            
            hist.at(i)->SetStats(kFALSE);
            hist.at(i)     ->SetMarkerStyle(20);
            hist.at(i)     ->SetMarkerSize(0.5);
        } 
        
        // Scale EXT
        else if (i == _util.k_ext){
            
            hist.at(i)->SetStats(kFALSE);
            // hist.at(i)->Scale(intime_scale_factor);

        }

        // Scale Dirt
        else if (i == _util.k_dirt){

            hist.at(i)->SetStats(kFALSE);
            // hist.at(i)->Scale(dirt_scale_factor);
            hist.at(i)     ->SetFillColor(2);
            hist.at(i)     ->SetFillStyle(3354);
            hist_integrals.at(i) = hist.at(i)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_dirt);

        }
        
        // Scale MC
        else {
            hist.at(i)->SetStats(kFALSE);
            // hist.at(i)->Scale(mc_scale_factor);
            hist.at(i)->SetFillColor(30);
            hist_integrals.at(i) = hist.at(i)->Integral();
            integral_mc_ext += hist_integrals.at(i);
        }

    }

    // Subtract the off beam from the on beam
    hist.at(_util.k_data)->Add(hist.at(_util.k_ext), -1);
    hist_integrals.at(_util.k_data) = hist.at(_util.k_data)->Integral();

    // Add the histograms to the stack
    h_stack->Add(hist.at(_util.k_mc));
    h_stack->Add(hist.at(_util.k_dirt)); 

    
    h_stack->Draw("hist");
    hist.at(_util.k_data)->Draw("same PE");

    // h_stack->GetYaxis()->SetTitle("Entries");

    // MC error histogram ------------------------------------------------------
    TH1D * h_error_hist = (TH1D*) hist.at(_util.k_mc)->Clone("h_error_hist");    
    h_error_hist->Add(hist.at(_util.k_dirt), 1);
    
    
    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2 hist same");    

    TLegend *leg_stack = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(hist.at(_util.k_data), "Data - EXT",   "lep");
    leg_stack->AddEntry(hist.at(_util.k_dirt), "Dirt",   "f");
    leg_stack->AddEntry(hist.at(_util.k_mc),   "Overlay","f");


    leg_stack->Draw();

    bottomPad->cd();

    h_ratio      = (TH1D*) hist.at(_util.k_data)->Clone("h_ratio");
    h_mc_ext_sum = (TH1D*) hist.at(_util.k_mc)  ->Clone("h_mc_ext_sum");
    
    // Add the dirt to overlay
    h_mc_ext_sum->Add(hist.at(_util.k_dirt), 1);
    

    h_ratio->Add(h_mc_ext_sum, -1);
    h_ratio->Divide(h_mc_ext_sum);

    h_ratio->GetXaxis()->SetLabelSize(12);
    h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetYaxis()->SetLabelSize(11);
    h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetXaxis()->SetTitleOffset(3.0);
    h_ratio->GetXaxis()->SetTitleSize(17);
    h_ratio->GetXaxis()->SetTitleFont(46);
    h_ratio->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

    h_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);
    h_ratio->GetYaxis()->SetTitle("(Data - MC) / MC ");
    h_ratio->GetYaxis()->SetTitleSize(13);
    h_ratio->GetYaxis()->SetTitleFont(44);
    h_ratio->GetYaxis()->SetTitleOffset(1.5);
    h_ratio->SetTitle(" ");
    h_ratio->Draw("E");

    // Draw the run period on the plot
    Draw_Run_Period(c);

    // Draw Data to MC ratio
    Draw_Data_MC_Ratio(c, double(hist_integrals.at(_util.k_data)*1.0/integral_mc_ext*1.0) );
    
    // Draw other data specifc quantities
    Draw_Data_POT(c, Data_POT);

    // Add the weight labels
    Draw_WeightLabels(c);
    
    c->Print(print_name);

}
// -----------------------------------------------------------------------------
void histogram_plotter::MakeEfficiencyPlot(const char* print_name){

    TTree *mc_tree;

    TFile *f_mc;

    // The uglyest hardcoded monstosity known to the universe. SORT this out KRISH...
    f_mc = TFile::Open("files/nuexsec_mc_run1.root");

    std::vector<double> *efficiency_v = nullptr; // efficiency vector
    std::vector<double> *purity_v     = nullptr; // purity vector

    _util.GetTree(f_mc, mc_tree, "mc_tree_out");
    mc_tree->SetBranchAddress("efficiency_v", &efficiency_v);
    mc_tree->SetBranchAddress("purity_v",     &purity_v);

    mc_tree->GetEntry(0); 

    TCanvas *c = new TCanvas();
    TH1D* h_eff = new TH1D("h_efficiency", "", efficiency_v->size(), 0, efficiency_v->size());
    TH1D* h_pur = new TH1D("h_purity", "", efficiency_v->size(), 0, efficiency_v->size());

    // c->SetGrid();
    c->SetGridy();

    TLegend *leg_stack = new TLegend(0.8, 0.91, 0.90, 0.72);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    
    
    
    for (unsigned int k=0; k < efficiency_v->size();k++){
        h_eff ->Fill(_util.cut_dirs.at(k).c_str(), efficiency_v->at(k));
        h_pur ->Fill(_util.cut_dirs.at(k).c_str(), purity_v->at(k));
        h_eff->SetBinError(k+1, 0);
        h_pur->SetBinError(k+1, 0);
    }
    
    leg_stack->AddEntry(h_eff, "Efficiency","lp");
    leg_stack->AddEntry(h_pur, "Purity",    "lp");

    h_eff->GetYaxis()->SetRangeUser(0,1.1);
    h_eff->SetStats(kFALSE);
    h_eff->SetMarkerStyle(20);
    h_eff->SetMarkerSize(0.5);
    h_eff->SetLineWidth(2);

    h_eff->Draw("LP");

    h_pur->SetLineColor(kRed+2);
    h_pur->SetStats(kFALSE);
    h_pur->SetMarkerStyle(20);
    h_pur->SetMarkerSize(0.5);
    h_pur->SetLineWidth(2);
    h_pur->Draw("LP,same");

    leg_stack->Draw();

    // Draw vertical lines to help the eye
    TLine *line;
    for (unsigned int l=1; l < efficiency_v->size()+1; l++){
        line  = new TLine( h_eff->GetBinCenter(l) ,   0 , h_eff->GetBinCenter(l)  ,  1.1);
        line->SetLineColor(12);
        line->SetLineStyle(kDotted);
        line->Draw();
    }

    // Draw the run period on the plot
    Draw_Run_Period(c);

    c->Print(print_name);


}
// -----------------------------------------------------------------------------
void histogram_plotter::MakeEfficiencyPlotByCut(const char* print_name,const char *run_period ){

    std::vector<TH1D*>  hist(_util.k_cuts_MAX); // The vector of histograms from the file for the plot
    std::vector<TEfficiency*> TEff_v(_util.k_cuts_MAX);

    TH1D* h_clone;

    TCanvas * c;

    // Loop over the classifications and get the histograms
    for (unsigned int i=0; i <_util.k_cuts_MAX; i++){

        // MC
        _util.GetHist(f_nuexsec, hist.at(i), Form("TEff/h_true_nu_E_%s", _util.cut_dirs.at(i).c_str()));
        if (hist.at(i) == NULL) return;
    }

    
    
    for (int p = 1; p < _util.k_cuts_MAX; p++){

        c = new TCanvas();

        // TEff_v.at(p) = new TEfficiency(*hist.at(p), *hist.at(_util.k_unselected));
        // TEff_v.at(p)->Draw("AP,same");
        h_clone = (TH1D*)hist.at(p)->Clone("h_clone");
        h_clone->Sumw2();
        h_clone->Divide(hist.at(_util.k_unselected));

        h_clone->SetStats(kFALSE);
        h_clone->SetTitle(Form("%s;True #nu_{e} Energy [GeV]; Efficiency", _util.cut_dirs.at(p).c_str()));
        h_clone->GetYaxis()->SetRangeUser(0,1);
        h_clone->SetLineColor(kBlack);
        h_clone->SetLineWidth(2);
        h_clone->Draw("E same");

        TH1D* h_true_nue = (TH1D*)hist.at(_util.k_unselected)->Clone("h_clone_true_nue");

        Float_t rightmax = 1.1*h_true_nue->GetMaximum();
        Float_t scale = gPad->GetUymax()/rightmax;
        h_true_nue ->SetLineColor(kRed+2);
        h_true_nue ->SetLineWidth(2);
        h_true_nue->Scale(scale);
        h_true_nue->Draw("hist,same");

        TGaxis *axis = new TGaxis(gPad->GetUxmax()+3,gPad->GetUymin(),gPad->GetUxmax()+3, gPad->GetUymax(),0,rightmax,510,"+L");
        axis->SetLineColor(kRed+2);
        axis->SetLabelColor(kRed+2);
        axis->SetTitle("True #nu_{e} Events in FV");
        axis->Draw();


        c->Print(Form("plots/run%s/Efficiency/TEff_%s.pdf", run_period, _util.cut_dirs.at(p).c_str() ));

    }
    
   



}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
