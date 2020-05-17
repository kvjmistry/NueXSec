#include "../include/SystematicsHelper.h"

// -----------------------------------------------------------------------------
void SystematicsHelper::Initialise(const char *_run_period, utility _utility){

    std::cout << "Initalising Systematics Helper..." << std::endl;
    _util = _utility;

    // To be added in as a configurable parameter
    std::string _mode = "ext";

    // Off beam mode to compare bnb and numi off beam samples
    if (_mode == "ext"){
        var_string = { "BNB", "NuMI" };
        mode = _mode;
    }

    // Get the POT of the variations from the file
    GetPOT(_run_period);

    // Set the run period
    run_period = std::string(_run_period);

    // Resize the file vector
    f_vars.resize(k_vars_MAX);

    // Get the variation files
    for (unsigned int l =0; l < var_string.size(); l++){
        
        // Standard variation mode
        if (mode == "default")  {
            f_vars.at(l) = new TFile( Form("files/nuexsec_run%s_%s_merged.root", _run_period, var_string.at(l).c_str() ), "READ");
        }
        // Off beam mode
        else if (mode == "ext") {
            f_vars.at(l) = new TFile( Form("files/nuexsec_ext_run%s_%s.root", _run_period, var_string.at(l).c_str() ), "READ");
        }
        else {
            std::cout << "Error I dont know what mode you have configured..." << mode << std::endl;
            exit(1);
        }
        
        
        // Maybe add something here to pickup non processed variation
    }

    // Now loop over events and caluclate the cross section
    MakeHistograms();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::GetPOT(const char* run_period){

    std::cout << "Getting the POT for the variation files" << std::endl;

    POT_v.resize(var_string.size(), 1.0);

    std::string line;

    std::string varname;
    std::string value;

    std::string pot_mode = "_MC_POT_"; // default mode

    if (mode == "ext") pot_mode =  "_EXT_trig_";

    std::string POT_run_config = "Run" + std::string(run_period) + pot_mode;
    
    // Loop over the config ist
    for (unsigned int p = 0; p < var_string.size(); p++){

        std::ifstream myfile ("config.txt");

        if (myfile.is_open()) {
            
            // std::cout << var_string.at(p) <<  std::endl;

            // Loop over lines in file
            while ( getline (myfile,line) ) {

                std::istringstream ss(line);
                ss >> varname >> value;

                // Found the correct variation file 
                std::string POT_name_match = POT_run_config + var_string.at(p);

                if (varname == POT_name_match ) {
                    std::cout << "Found match for: " << varname << " "<< std::stod(value) <<  std::endl;
                    POT_v.at(p)= std::stod(value);
                    break;
                }
                
            }
           
        }
        else std::cout << "Unable to open config file, bad things are going to happen..." << std::endl; 

        myfile.close();
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetVariationProperties(TH1D* h, int index){

    if (mode == "default"){
        if (index == k_CV){
            h->SetLineColor(kBlack);
        }
        else if (index == k_bnb_diffusion){
            h->SetLineColor(kAzure-5);
        }
    }
    else {
        
        if (index == k_NuMI){
            h->SetLineColor(kBlack);
        }
        else if (index == k_BNB){
            h->SetLineColor(kAzure-5);
        }

    }
    
    h->SetLineWidth(2);
    h->GetYaxis()->SetTitle("Entries");

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetTPadOptions(TPad * topPad, TPad * bottomPad ){

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
void SystematicsHelper::CreateDirectory(std::string folder, std::string run_period){

    std::string a = "if [ ! -d \"plots/";
    std::string b = "run" + std::string(run_period) + "/" + folder;
    std::string c = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
    std::string d = "run" + std::string(run_period) + "/" + folder;
    std::string e = "; fi";
    std::string command = a + b + c + d + e ;
    system(command.c_str()); 

}
// -----------------------------------------------------------------------------
void SystematicsHelper::MakeHistograms(){
    
    for (unsigned int i = 0 ; i < _util.k_cuts_MAX; i++){
        
        // Default detector systematics mode
        if (mode == "defaut"){
            // Create the directory
            CreateDirectory("/detvar/comparisons/cuts/" + _util.cut_dirs.at(i), run_period);

            // Space Charge Corrected X position comparision plot
            PlotVariations("h_reco_vtx_x_sce", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_vtx_x_sce.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex X [cm]");

            // Space Charge Corrected Y position comparision plot
            PlotVariations("h_reco_vtx_y_sce", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_vtx_y_sce.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Y [cm]");


            // Space Charge Corrected X position comparision plot
            PlotVariations("h_reco_vtx_z_sce", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_vtx_z_sce.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Z [cm]");
        
            // Flash Time
            PlotVariations("h_reco_flash_time", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_flash_time.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Flash Time [#mus]");

            // Leading Shower Phi
            PlotVariations("h_reco_leading_shower_phi", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_leading_shower_phi.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Phi [degrees]");

            // Leading Shower Theta
            PlotVariations("h_reco_leading_shower_theta", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_leading_shower_theta.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Theta [degrees]");


            // Shower Multiplicty
            PlotVariations("h_reco_shower_multiplicity", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_shower_multiplicity.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Multiplicty");

            // Track Multiplicty
            PlotVariations("h_reco_track_multiplicity", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_track_multiplicity.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Track Multiplicity");

            
            // Topological Score
            PlotVariations("h_reco_topological_score", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_topological_score.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Topological Score");

            
            // Shower Track Fitter dedx Y plane
            PlotVariations("h_reco_shr_tkfit_dedx_y", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_shr_tkfit_dedx_y.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Collection Plane dEdx (track fitter) [MeV/cm]");

            // Reco Electron Neutrino E
            PlotVariations("h_reco_nu_e", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_nu_e.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Neutrino Energy [GeV]");

            // Leading Shower Energy
            PlotVariations("h_reco_shower_energy_tot_cali", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_shower_energy_tot_cali.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Leading Shower Energy [GeV]");

        }
        // Ext mode
        else if (mode == "ext"){

            // Create the directory
            CreateDirectory("/ext/comparisons/cuts/" + _util.cut_dirs.at(i), run_period);

            // Space Charge Corrected X position comparision plot
            PlotVariationsEXT("h_reco_vtx_x_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_x_sce.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex X [cm]");

        }
        // dont know what your trying to configure tbh ;)
        else {
            std::cout << "something is going wrong here..." << std::endl;
        }
        
    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotVariations(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name ){

    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    hist.resize(k_vars_MAX);
    hist_ratio.resize(k_vars_MAX);

    TCanvas * c      = new TCanvas();
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    SetTPadOptions(topPad, bottomPad );

    // Loop over the variations and get the histograms
    for (unsigned int k=0; k < f_vars.size(); k++){
        
        // Loop over the classifications and get the histograms
        for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

            // Only want the MC piece -- may want to add in dirt too? -- will need to separately scale that hitogram though
            if ( i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt ) continue;

            // Get all the MC histograms and add them together
            TH1D *h_temp;

            _util.GetHist(f_vars.at(k), h_temp, Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            
            // First case so clone the histogram
            if (i == 0) hist.at(k) = (TH1D*) h_temp->Clone("h_sum_hist");
            else hist.at(k)->Add(h_temp, 1);
        }
        
    }

    // Now scale the histograms to POT
    for (unsigned int y=0; y < hist.size(); y++ ){
        double scale_fact = POT_v.at(k_CV) / POT_v.at(y);
        // std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        if (y == k_CV){
            // Clone a histogram to plot the CV error as a grey band
            h_error_hist = (TH1D*) hist.at(k_CV)->Clone("h_error_hist");
            h_error_hist->SetFillColorAlpha(12, 0.15);
            h_error_hist->SetLineWidth(2);
            h_error_hist->SetLineColor(kBlack);
            h_error_hist->GetYaxis()->SetTitle("Entries");
            h_error_hist->GetYaxis()->SetTitleFont(46);
            h_error_hist->GetYaxis()->SetTitleSize(13);
            h_error_hist->Draw("E2");
        }


        // Save clones of the histograms for doing the ratios
        hist_ratio.at(y) = (TH1D*) hist.at(y)->Clone(Form("h_ratio_%s", var_string.at(y).c_str()));
        hist_ratio.at(y)->Divide(hist.at(k_CV)); 

        // Set the customisation of the histogram
        SetVariationProperties(hist.at(y), y);

        // Draw the histograms
        if (y == k_CV) hist.at(y)->Draw("hist, same");
        else hist.at(y)->Draw("hist,E, same");
    }

    // Legend
    TLegend *leg = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_error_hist, "CV",   "lf");
    leg->AddEntry(hist.at(k_bnb_diffusion), "BNB Diffusion", "l");
    leg->Draw();

    // Now draw the ratio
    bottomPad->cd();

    for (unsigned int k =0; k < hist_ratio.size(); k++){

        SetVariationProperties(hist_ratio.at(k), k);
        hist_ratio.at(k)->SetLineWidth(1.5);


        hist_ratio.at(k)->GetXaxis()->SetLabelSize(12);
        hist_ratio.at(k)->GetXaxis()->SetLabelFont(43); 
        hist_ratio.at(k)->GetYaxis()->SetLabelSize(11);
        hist_ratio.at(k)->GetYaxis()->SetLabelFont(43);
        hist_ratio.at(k)->GetXaxis()->SetTitleOffset(3.0);
        hist_ratio.at(k)->GetXaxis()->SetTitleSize(17);
        hist_ratio.at(k)->GetXaxis()->SetTitleFont(46);
        hist_ratio.at(k)->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        hist_ratio.at(k)->GetYaxis()->SetRangeUser(0., 2.0);
        hist_ratio.at(k)->GetYaxis()->SetTitle("Variation / CV");
        hist_ratio.at(k)->GetYaxis()->SetTitleSize(13);
        hist_ratio.at(k)->GetYaxis()->SetTitleFont(44);
        hist_ratio.at(k)->GetYaxis()->SetTitleOffset(1.5);
        hist_ratio.at(k)->SetTitle(" ");
        hist_ratio.at(k)->GetXaxis()->SetTitle(x_axis_name);
        
        if (k == 0) hist_ratio.at(k)->Draw("hist,same");
        else hist_ratio.at(k)->Draw("hist,E,same");
    }

    // Draw the error hist 
    TH1D* h_ratio_error = (TH1D*) h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");


    // Draw the run period on the plot
    // Draw_Run_Period(c);

    // Add the weight labels
    // Draw_WeightLabels(c);
    
    c->Print(print_name);



}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotVariationsEXT(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name ){

    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    
    hist.resize(var_string.size());
    hist_ratio.resize(var_string.size());

    TCanvas * c      = new TCanvas();
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    SetTPadOptions(topPad, bottomPad );


    // Loop over the variations and get the histograms
    for (unsigned int k=0; k < f_vars.size(); k++){
        
        // Loop over the classifications and get the histograms
        for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

            // Only want the EXT piece
            if ( i != _util.k_leg_ext ) continue;

            _util.GetHist(f_vars.at(k), hist.at(k), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            
        }
        
    }


    // Now scale the histograms to POT
    for (unsigned int y=0; y < hist.size(); y++ ){
        double scale_fact = POT_v.at(k_NuMI) / POT_v.at(y);
        std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        // Do area norm
        hist.at(y)->Scale(hist.at(k_NuMI)->Integral() / hist.at(y)->Integral() );

        if (y == k_NuMI){
            // Clone a histogram to plot the CV error as a grey band
            h_error_hist = (TH1D*) hist.at(k_NuMI)->Clone("h_error_hist");
            h_error_hist->SetFillColorAlpha(12, 0.15);
            h_error_hist->SetLineWidth(2);
            h_error_hist->SetLineColor(kBlack);
            h_error_hist->GetYaxis()->SetTitle("Entries");
            h_error_hist->GetYaxis()->SetTitleFont(46);
            h_error_hist->GetYaxis()->SetTitleSize(13);
            h_error_hist->Draw("E2");
        }


        // Save clones of the histograms for doing the ratios
        hist_ratio.at(y) = (TH1D*) hist.at(y)->Clone(Form("h_ratio_%s", var_string.at(y).c_str()));
        hist_ratio.at(y)->Divide(hist.at(k_NuMI)); 

        // Set the customisation of the histogram
        SetVariationProperties(hist.at(y), y);

        // Draw the histograms
        if (y == k_NuMI) hist.at(y)->Draw("hist, same");
        else hist.at(y)->Draw("hist,E, same");

    }

    // Legend
    TLegend *leg = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_error_hist, "NuMI EXT",   "lf");
    leg->AddEntry(hist.at(k_BNB), "BNB EXT", "l");
    leg->Draw();

    // Now draw the ratio
    bottomPad->cd();

    for (unsigned int k =0; k < hist_ratio.size(); k++){

        SetVariationProperties(hist_ratio.at(k), k);
        hist_ratio.at(k)->SetLineWidth(1.5);


        hist_ratio.at(k)->GetXaxis()->SetLabelSize(12);
        hist_ratio.at(k)->GetXaxis()->SetLabelFont(43); 
        hist_ratio.at(k)->GetYaxis()->SetLabelSize(11);
        hist_ratio.at(k)->GetYaxis()->SetLabelFont(43);
        hist_ratio.at(k)->GetXaxis()->SetTitleOffset(3.0);
        hist_ratio.at(k)->GetXaxis()->SetTitleSize(17);
        hist_ratio.at(k)->GetXaxis()->SetTitleFont(46);
        hist_ratio.at(k)->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        hist_ratio.at(k)->GetYaxis()->SetRangeUser(0., 2.0);
        hist_ratio.at(k)->GetYaxis()->SetTitle("Variation / CV");
        hist_ratio.at(k)->GetYaxis()->SetTitleSize(13);
        hist_ratio.at(k)->GetYaxis()->SetTitleFont(44);
        hist_ratio.at(k)->GetYaxis()->SetTitleOffset(1.5);
        hist_ratio.at(k)->SetTitle(" ");
        hist_ratio.at(k)->GetXaxis()->SetTitle(x_axis_name);
        
        if (k == 0) hist_ratio.at(k)->Draw("hist,same");
        else hist_ratio.at(k)->Draw("hist,E,same");
    }

    // Draw the error hist 
    TH1D* h_ratio_error = (TH1D*) h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");


    // Draw the run period on the plot
    // Draw_Run_Period(c);

    // Add the weight labels
    // Draw_WeightLabels(c);
    
    c->Print(print_name);



}
// -----------------------------------------------------------------------------