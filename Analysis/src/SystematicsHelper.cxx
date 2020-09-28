#include "../include/SystematicsHelper.h"

// -----------------------------------------------------------------------------
void SystematicsHelper::Initialise(Utility _utility){

    std::cout << "Initalising Systematics Helper..." << std::endl;
    _util = _utility;

    // Open the outoput systematics file
    file_sys_var = TFile::Open(Form("files/run%s_sys_var.root", _util.run_period),"UPDATE");

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
    
    // Get the POT of the variations from the file
    GetPOT();

    // Off beam mode to compare bnb and numi off beam samples
    if (std::string(_util.sysmode) == "ext"){
        var_string = { "NuMI", "BNB" };
    }

    // If we choose this mode then we actually want to use a different initialiser
    if (std::string(_util.sysmode) == "reweight"){
        InitialiseReweightingMode();
        return;
    }

    // We want to get the systematic uncertainty for specific plots in each cut
    if (std::string(_util.sysmode) == "reweightcuts"){
        InitialiseReweightingModeCut();
        return;
    }

    // Resize the file vector
    f_vars.resize(k_vars_MAX);

    // Get the variation files
    for (unsigned int l =0; l < var_string.size(); l++){
        
        // Standard variation mode
        if (std::string(_util.sysmode) == "default")  {
            f_vars.at(l) = new TFile( Form("files/nuexsec_run%s_%s_merged.root", _util.run_period, var_string.at(l).c_str() ), "READ");
        }
        // Off beam mode
        else if (std::string(_util.sysmode) == "ext") {
            f_vars.at(l) = new TFile( Form("files/nuexsec_ext_run%s_%s.root", _util.run_period, var_string.at(l).c_str() ), "READ");
        }
        else {
            std::cout << "Error I dont know what mode you have configured..." << std::string(_util.sysmode) << std::endl;
            exit(1);
        }
        
        
        // Maybe add something here to pickup non processed variation
    }

    // Now loop over events and caluclate the cross section
    MakeHistograms();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::GetPOT(){

    std::cout << "Getting the POT for the variation files" << std::endl;

    POT_v.resize(var_string.size(), 1.0);

    std::string line;

    std::string varname;
    std::string value;

    std::string pot_mode = "_MC_POT_"; // default mode

    if (mode == "ext") pot_mode =  "_EXT_trig_";

    std::string POT_run_config = "Run" + std::string(_util.run_period) + pot_mode;
    
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
void SystematicsHelper::Draw_Area_Norm(TCanvas* c){
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
void SystematicsHelper::SetRatioOptions(TH1D* hist ){

    hist->GetXaxis()->SetLabelSize(0.13);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetTitleSize(0.13);
    hist->GetYaxis()->SetLabelSize(0.13);
    hist->GetYaxis()->SetRangeUser(-50, 50);
    hist->GetYaxis()->SetTitleSize(12);
    hist->GetYaxis()->SetTitleFont(44);
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleOffset(1.5);
    hist->SetTitle(" ");

}
// -----------------------------------------------------------------------------
void SystematicsHelper::MakeHistograms(){
    
    for (unsigned int i = 0 ; i < _util.k_cuts_MAX; i++){
        
        // Default detector systematics mode
        if (mode == "default"){
  
            // Create the directory
            _util.CreateDirectory("/detvar/comparisons/cuts/" + _util.cut_dirs.at(i));

            // Create the directory for sysvar
            _util.CreateDirectory("/systvar/comparisons/cuts/" + _util.cut_dirs.at(i));

            for(unsigned int j=0; j < _util.vec_hist_name.size(); j++){

                SysVariations(Form("%s", _util.vec_hist_name.at(j).c_str()), Form("plots/run%s/systvar/comparisons/cuts/%s/%s.pdf", _util.run_period, _util.cut_dirs.at(i).c_str(), _util.vec_hist_name.at(j).c_str()),
                            _util.cut_dirs.at(i), _util.vec_axis_label.at(j).c_str(), _util.cut_dirs.at(i).c_str(), _util.vec_hist_name.at(j).c_str());

                PlotVariations(Form("%s", _util.vec_hist_name.at(j).c_str()), Form("plots/run%s/detvar/comparisons/cuts/%s/%s.pdf", _util.run_period, _util.cut_dirs.at(i).c_str(), _util.vec_hist_name.at(j).c_str()),
                            _util.cut_dirs.at(i), _util.vec_axis_label.at(j).c_str());
            }
        }
        // Ext mode
        else if (mode == "ext"){

            // Create the directory
            _util.CreateDirectory("/ext/comparisons/cuts/" + _util.cut_dirs.at(i));

            // Space Charge Corrected X position comparision plot
            PlotVariationsEXT("h_reco_vtx_x_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_x_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex X [cm]");

            // Space Charge Corrected Y position comparision plot
            PlotVariationsEXT("h_reco_vtx_y_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_y_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Y [cm]");


            // Space Charge Corrected X position comparision plot
            PlotVariationsEXT("h_reco_vtx_z_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_z_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Z [cm]");
        
            // Leading Shower Phi
            PlotVariationsEXT("h_reco_leading_shower_phi", Form("plots/run%s/ext/comparisons/cuts/%s/reco_leading_shower_phi.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Phi [degrees]");

            // Leading Shower Theta
            PlotVariationsEXT("h_reco_leading_shower_theta", Form("plots/run%s/ext/comparisons/cuts/%s/reco_leading_shower_theta.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Theta [degrees]");


            // Shower Multiplicty
            PlotVariationsEXT("h_reco_shower_multiplicity", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shower_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Multiplicty");

            // Track Multiplicty
            PlotVariationsEXT("h_reco_track_multiplicity", Form("plots/run%s/ext/comparisons/cuts/%s/reco_track_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Track Multiplicity");

            
            // Topological Score
            PlotVariationsEXT("h_reco_topological_score", Form("plots/run%s/ext/comparisons/cuts/%s/reco_topological_score.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Topological Score");

            
            // Shower Track Fitter dedx Y plane
            PlotVariationsEXT("h_reco_shr_tkfit_dedx_y", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shr_tkfit_dedx_y.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Collection Plane dEdx (track fitter) [MeV/cm]");

            // Reco Electron Neutrino E
            PlotVariationsEXT("h_reco_nu_e", Form("plots/run%s/ext/comparisons/cuts/%s/reco_nu_e.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Neutrino Energy [GeV]");

            // Leading Shower Energy
            PlotVariationsEXT("h_reco_shower_energy_tot_cali", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shower_energy_tot_cali.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Leading Shower Energy [GeV]");

            // Flash PE
            PlotVariationsEXT("h_reco_flash_pe", Form("plots/run%s/ext/comparisons/cuts/%s/reco_flash_pe.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Flash PE [PE]");

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

    _util.SetTPadOptions(topPad, bottomPad );

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
        hist_ratio.at(k)->SetLineWidth(1);


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
    // _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);
    
    c->Print(print_name);



}
// ----------------------------------------------------------------------------
void SystematicsHelper::SysVariations(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name, std::string folder_name, std::string plot_name){

    // last updated on Sept 18 by Marina Reggiani-Guzzo

    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_diff; // The vector of histogram differentes between CV and the vatiation (variation-CV)
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    hist.resize(k_vars_MAX);
    hist_diff.resize(k_vars_MAX);
    hist_ratio.resize(k_vars_MAX);

    TCanvas * c      = new TCanvas();
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    _util.SetTPadOptions(topPad, bottomPad );

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

        // calculate difference between cv and the histogram with variation
        hist_diff.at(y) = (TH1D*) hist.at(y)->Clone(Form("h_diff_%s", var_string.at(y).c_str()));
        hist_diff.at(y)->Add(hist.at(k_CV),-1); 

        // Save clones of the histograms for doing the ratios
        hist_ratio.at(y) = (TH1D*) hist_diff.at(y)->Clone(Form("h_ratio_%s", var_string.at(y).c_str()));
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
        hist_ratio.at(k)->SetLineWidth(1);


        hist_ratio.at(k)->GetXaxis()->SetLabelSize(12);
        hist_ratio.at(k)->GetXaxis()->SetLabelFont(43); 
        hist_ratio.at(k)->GetYaxis()->SetLabelSize(11);
        hist_ratio.at(k)->GetYaxis()->SetLabelFont(43);
        hist_ratio.at(k)->GetXaxis()->SetTitleOffset(3.0);
        hist_ratio.at(k)->GetXaxis()->SetTitleSize(17);
        hist_ratio.at(k)->GetXaxis()->SetTitleFont(46);
        hist_ratio.at(k)->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        hist_ratio.at(k)->GetYaxis()->SetRangeUser(-3.0, 3.0);
        hist_ratio.at(k)->GetYaxis()->SetTitle("(Variation-CV) / CV");
        hist_ratio.at(k)->GetYaxis()->SetTitleSize(13);
        hist_ratio.at(k)->GetYaxis()->SetTitleFont(44);
        hist_ratio.at(k)->GetYaxis()->SetTitleOffset(1.5);
        hist_ratio.at(k)->SetTitle(" ");
        hist_ratio.at(k)->GetXaxis()->SetTitle(x_axis_name);
        
        if (k == 0) hist_ratio.at(k)->Draw("hist,same");
        else {

        hist_ratio.at(k)->Draw("hist,E,same");

        }

    }

    // -----------------------------------------------------------------
    // calculate and save the total detector systematics uncertainty 

    // create a temporary histogram that will be used to calculate the total detector sys
    TH1D *h_det_sys_tot;

    file_sys_var->cd();
    
    // loop over variations for a given variable
    for (unsigned int k = 0; k < f_vars.size(); k++){
        
        // ---- save the histograms into different directories inside the root file

        if(!file_sys_var->GetDirectory(Form("%s/%s", folder_name.c_str(), var_string.at(k).c_str()))) {
            file_sys_var->mkdir(Form("%s/%s", folder_name.c_str(), var_string.at(k).c_str())); // if the directory does not exist, create it
        }

        file_sys_var->cd(Form("%s/%s", folder_name.c_str(), var_string.at(k).c_str())); // open the directory
    
        hist_ratio.at(k)->SetDirectory(gDirectory); // set in which dir the hist_ratio.at(k) is going to be written
        hist_ratio.at(k)->Write(Form("%s", plot_name.c_str()), TObject::kOverwrite);  // write the histogram to the file
    
        // ---- on the same go, calculate the total detector sys uncertainty

        if( k == 0 ){
            // this is the first histogram, just square it and write it to file
            h_det_sys_tot = (TH1D*) hist_ratio.at(k)->Clone();
            h_det_sys_tot->Multiply(h_det_sys_tot); // square the histogram
        }

        else{
            // this is not the first histogram
            hist_ratio.at(k)->Multiply(hist_ratio.at(k)); // first square the histogram
            h_det_sys_tot->Add(hist_ratio.at(k)); // add it to the existing histogram
        }

    }

    // calculate the square root of the histogram before saving it to the file
    for(int k = 1; k <= h_det_sys_tot->GetNbinsX(); k++){
        h_det_sys_tot->SetBinContent( k , TMath::Sqrt(h_det_sys_tot->GetBinContent(k)) );
    }

    // save the total detector sys uncertainty in the file
    
    if(!file_sys_var->GetDirectory(Form("%s/TotalDetectorSys", folder_name.c_str()))) {
                file_sys_var->mkdir(Form("%s/TotalDetectorSys", folder_name.c_str())); // if the directory does not exist, create it
    }

    file_sys_var->cd(Form("%s/TotalDetectorSys", folder_name.c_str())); // open the directory
    h_det_sys_tot->SetDirectory(gDirectory);
    h_det_sys_tot->Write(Form("%s", plot_name.c_str()), TObject::kOverwrite); 
  
    // -----------------------------------------------------------------


  /*
    TH1D* h_ratio_error = (TH1D*) h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");
    */


    // Draw the run period on the plot
    // _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);
    
    c->Print(print_name);

}
// ----------------------------------------------------------------------------
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

    _util.SetTPadOptions(topPad, bottomPad );


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
        // std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        // Do area norm
        std::cout << "BNB to NuMI scale Factor: " << hist.at(k_NuMI)->Integral() / hist.at(y)->Integral() << std::endl;
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
        hist_ratio.at(k)->SetLineWidth(1);


        hist_ratio.at(k)->GetXaxis()->SetLabelSize(12);
        hist_ratio.at(k)->GetXaxis()->SetLabelFont(43); 
        hist_ratio.at(k)->GetYaxis()->SetLabelSize(11);
        hist_ratio.at(k)->GetYaxis()->SetLabelFont(43);
        hist_ratio.at(k)->GetXaxis()->SetTitleOffset(3.0);
        hist_ratio.at(k)->GetXaxis()->SetTitleSize(17);
        hist_ratio.at(k)->GetXaxis()->SetTitleFont(46);
        hist_ratio.at(k)->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        hist_ratio.at(k)->GetYaxis()->SetRangeUser(0., 2.0);
        hist_ratio.at(k)->GetYaxis()->SetTitle("NuMI / BNB");
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
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Draw area normalisation
    Draw_Area_Norm(c);
    
    c->Print(print_name);



}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseReweightingMode(){

    gStyle->SetOptStat(0);

    // Load in the input file
    // Should we add more protection to this command??
    f_nuexsec = new TFile( Form("files/crosssec_run%s.root", _util.run_period ), "READ");

    InitialsePlotCV();

    // Now lets initialse the vectors which will store the total uncertainties
    InitialiseUncertaintyVectors();

    // Initialise the covariance matrices
    int n_bins = cv_hist_vec.at(k_var_reco_el_E).at(0)->GetNbinsX();
    h_cov_tot         = new TH2D("h_cov_tot",         "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_sys         = new TH2D("h_cov_sys",         "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_stat        = new TH2D("h_cov_stat",        "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_genie_uni   = new TH2D("h_cov_genie_uni",   "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_genie_multi = new TH2D("h_cov_genie_multi", "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_hp          = new TH2D("h_cov_hp",          "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_beamline    = new TH2D("h_cov_beamline",    "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_dirt        = new TH2D("h_cov_dirt",        "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_pot         = new TH2D("h_cov_pot",         "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
    h_cov_reint       = new TH2D("h_cov_reint",       "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);


    // Loop over the cross-section variables
    for (unsigned int var = 0; var <  vars.size(); var++){

        // Comparison plots for data to MC
        CompareVariationXSec("RPA",              var, "RPA" );
        CompareVariationXSec("CCMEC",            var, "CC MEC" );
        CompareVariationXSec("AxFFCCQE",         var, "Ax FF CCQE" );
        CompareVariationXSec("VecFFCCQE",        var, "Vec FF CCQE" );
        CompareVariationXSec("DecayAngMEC",      var, "Decay Ang MEC" );
        CompareVariationXSec("ThetaDelta2Npi",   var, "Theta Delta 2N #pi" );
        CompareVariationXSec("ThetaDelta2NRad",  var, "Theta Delta 2N Rad" );
        CompareVariationXSec("RPA_CCQE_Reduced", var, "RPA CCQE Reduced" );
        CompareVariationXSec("NormCCCOH",        var, "Norm CC COH" );
        CompareVariationXSec("NormNCCOH",        var, "Norm NC COH" );

        // Plot the interaction unisims
        PlotReweightingModeUnisim("RPA",              var, "RPA" );
        PlotReweightingModeUnisim("CCMEC",            var, "CC MEC" );
        PlotReweightingModeUnisim("AxFFCCQE",         var, "Ax FF CCQE" );
        PlotReweightingModeUnisim("VecFFCCQE",        var, "Vec FF CCQE" );
        PlotReweightingModeUnisim("DecayAngMEC",      var, "Decay Ang MEC" );
        PlotReweightingModeUnisim("ThetaDelta2Npi",   var, "Theta Delta 2N #pi" );
        PlotReweightingModeUnisim("ThetaDelta2NRad",  var, "Theta Delta 2N Rad" );
        PlotReweightingModeUnisim("RPA_CCQE_Reduced", var, "RPA CCQE Reduced" );
        PlotReweightingModeUnisim("NormCCCOH",        var, "Norm CC COH" );
        PlotReweightingModeUnisim("NormNCCOH",        var, "Norm NC COH" );

        // Dirt
        PlotReweightingModeUnisim("Dirt",        var, "Dirt" );

        // POT Counting
        PlotReweightingModeUnisim("POT",        var, "POT Count." );

        // Plot the beamline unisims
        PlotReweightingModeUnisim("Horn1_x",            var, "Horn 1 x" );
        PlotReweightingModeUnisim("Horn_curr",          var, "Horn Current" );
        PlotReweightingModeUnisim("Horn1_y",            var, "Horn 1 y" );
        PlotReweightingModeUnisim("Beam_spot",          var, "Beam Spot Size" );
        PlotReweightingModeUnisim("Horn2_x",            var, "Horn 2 x" );
        PlotReweightingModeUnisim("Horn2_y",            var, "Horn 2 y" );
        PlotReweightingModeUnisim("Horn_Water",         var, "Horns Water" );
        PlotReweightingModeUnisim("Beam_shift_x",       var, "Beam shift x" );
        PlotReweightingModeUnisim("Beam_shift_y",       var, "Beam shift y" );
        PlotReweightingModeUnisim("Target_z",           var, "Target z" );
        PlotReweightingModeUnisim("Decay_pipe_Bfield",  var, "Decay pipe Bfield" );

        // Detector Variations
        PlotReweightingModeDetVar("LYRayleigh",                         var, k_LYRayleigh,                         var_string_pretty.at(k_LYRayleigh));
        PlotReweightingModeDetVar("SCE",                                var, k_SCE,                                var_string_pretty.at(k_SCE));
        PlotReweightingModeDetVar("LYAttenuation",                      var, k_LYAttenuation,                      var_string_pretty.at(k_LYAttenuation));
        PlotReweightingModeDetVar("Recomb2",                            var, k_Recomb2,                            var_string_pretty.at(k_Recomb2));
        PlotReweightingModeDetVar("WireModX",                           var, k_WireModX,                           var_string_pretty.at(k_WireModX));
        PlotReweightingModeDetVar("WireModYZ",                          var, k_WireModYZ,                          var_string_pretty.at(k_WireModYZ));
        PlotReweightingModeDetVar("WireModThetaXZ",                     var, k_WireModThetaXZ,                     var_string_pretty.at(k_WireModThetaXZ));
        PlotReweightingModeDetVar("WireModThetaYZ_withSigmaSplines",    var, k_WireModThetaYZ_withSigmaSplines,    var_string_pretty.at(k_WireModThetaYZ_withSigmaSplines));
        PlotReweightingModeDetVar("WireModThetaYZ_withoutSigmaSplines", var, k_WireModThetaYZ_withoutSigmaSplines, var_string_pretty.at(k_WireModThetaYZ_withoutSigmaSplines));
        PlotReweightingModeDetVar("WireModdEdX",                        var, k_WireModdEdX,                        var_string_pretty.at(k_WireModdEdX));

        // Plot the multisims
        PlotReweightingModeMultisim("weightsGenie", var,  "GENIE", 600);
        PlotReweightingModeMultisim("weightsReint", var,  "Geant Reinteractions", 1000);
        PlotReweightingModeMultisim("weightsPPFX",  var,  "PPFX", 600);
        PlotReweightingModeMultisim("MCStats",      var,  "MC Stats", 1000);
        
    }

    // Get the statistical uncertainties
    FillStatVector();

    // Compare the MC and Data X-Section
    CompareCVXSec();
    CompareCVXSecNoRatio();

    // Save the total covariance matrices
    _util.CreateDirectory("/Systematics/Covariance");
    SaveCovMatrix(h_cov_tot,         Form("plots/run%s/Systematics/Covariance/run%s_tot_cov.pdf",         _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_sys,         Form("plots/run%s/Systematics/Covariance/run%s_tot_sys_cov.pdf",     _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_stat,        Form("plots/run%s/Systematics/Covariance/run%s_tot_stat_cov.pdf",    _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_genie_uni,   Form("plots/run%s/Systematics/Covariance/run%s_genie_uni_cov.pdf",   _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_genie_multi, Form("plots/run%s/Systematics/Covariance/run%s_genie_multi_cov.pdf", _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_hp,          Form("plots/run%s/Systematics/Covariance/run%s_hp_cov.pdf",          _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_beamline,    Form("plots/run%s/Systematics/Covariance/run%s_beamline_cov.pdf",    _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_dirt,        Form("plots/run%s/Systematics/Covariance/run%s_dirt_cov.pdf",        _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_pot,         Form("plots/run%s/Systematics/Covariance/run%s_pot_cov.pdf",         _util.run_period, _util.run_period));
    SaveCovMatrix(h_cov_reint,       Form("plots/run%s/Systematics/Covariance/run%s_reint_cov.pdf",       _util.run_period, _util.run_period));

    // Create the directories
    _util.CreateDirectory("/Systematics/Beamline");
    _util.CreateDirectory("/Systematics/Genie_Unisim");
    
    // Plot the total beamline sys uncertainty
    PlotTotUnisim("Beamline");
    PlotTotUnisim("Genie_Unisim");

    // Print a summary of the results
    PrintUncertaintySummary();

    // Print the sqrt of the diagonals of the covariance matrix
    // loop over rows
    for (int row = 1; row < h_cov_sys->GetNbinsX()+1; row++) {
        
        // Loop over columns
        for (int col = 1; col < h_cov_sys->GetNbinsY()+1; col++) {
            
            // Only set the bin content of the diagonals
            double bin_diag = cv_hist_vec.at(k_var_reco_el_E).at(k_xsec_mcxsec)->GetBinContent(row);

            // 0.01 converts each percentage back to a number. We multiply this by the cv to get the deviate
            if (row == col) std::cout << 100 * std::sqrt(h_cov_sys->GetBinContent(row, col)) / bin_diag << std::endl;  
        }
    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetLabelName(std::string label, std::string &label_up, std::string &label_dn){

    if  (label == "Horn_curr"          ){
        label_up = "Horn_p2kA";
        label_dn = "Horn_m2kA";
    }
    else if  (label == "Horn1_x"       ){
        label_up = label + "_p3mm";
        label_dn = "Horm1_x_m3mm";
    }
    else if  (label == "Horn1_y"       ){
        label_up = label + "_p3mm";
        label_dn = label + "_m3mm";
    }
    else if  (label == "Beam_spot"    ){
        label_up = label + "_1_1mm";
        label_dn = label + "_1_5mm";
    }
    else if  (label == "Horn2_x"       ){
        label_up = label + "_p3mm";
        label_dn = "Horm2_x_m3mm";
    }
    else if  (label == "Horn2_y"       ){
        label_up = label + "_p3mm";
        label_dn = label + "_m3mm";
    }
    else if  (label == "Horn_Water"    ){
        label_up = "Horns_0mm_water";
        label_dn = "Horns_2mm_water";
    }
    else if  (label == "Beam_shift_x"  ){
        label_up = label + "_p1mm";
        label_dn = label + "_m1mm";
    }
    else if  (label == "Beam_shift_y"  ){
        label_up = label + "_p1mm";
        label_dn = label + "_m1mm";
    }
    else if  (label == "Target_z"      ){
        label_up = label + "_p7mm";
        label_dn = label + "_m7mm";
    }
    else if  (label == "Horn1_refined_descr" || label == "Decay_pipe_Bfield" || label == "Old_Horn_Geometry"){
        label_up = label;
        label_dn = label;
    }
    // Other variation
    else {
        label_up = label + "up";
        label_dn = label + "dn";
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeUnisim(std::string label, int var, std::string label_pretty){

    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe;
    
    // Resize to the number of universes
    h_universe.resize(2);
    h_universe.at(k_up).resize(xsec_types.size());
    h_universe.at(k_dn).resize(xsec_types.size());

    // Now get the histograms
    std::string label_up = label + "up";
    std::string label_dn = label + "dn";

    // Set the Unisim up down variation name
    SetLabelName(label, label_up, label_dn);

    // Check if its just a single on/off type variation
    // This case we dont want to plot the up/dn, but just once
    bool single_var = false;
    if (label_up == label_dn) single_var = true;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        _util.GetHist(f_nuexsec, h_universe.at(k_up).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), _util.run_period, label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, h_universe.at(k_dn).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), _util.run_period, label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
    }

    // Get the covariance matrices
    CalcMatrices(label, var, h_universe);

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c;
    
    // Now we want to draw them
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        c = new TCanvas("c", "c", 500, 500);
        topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        _util.SetTPadOptions(topPad, bottomPad);
        // topPad->SetRightMargin(0.10 );
        // bottomPad->SetRightMargin(0.10 );
        
        h_universe.at(k_up).at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));
        h_universe.at(k_up).at(k)->GetXaxis()->SetTitle("");
        h_universe.at(k_up).at(k)->GetXaxis()->SetLabelSize(0);

        h_universe.at(k_up).at(k)->Draw("hist");
        cv_hist_vec.at(var).at(k)->Draw("hist,same");
        if (!single_var) h_universe.at(k_dn).at(k)->Draw("hist,same");

        c->Update();

        double scale_val = h_universe.at(k_up).at(k)->GetMaximum();
        if (scale_val < cv_hist_vec.at(var).at(k)->GetMaximum()) scale_val = cv_hist_vec.at(var).at(k)->GetMaximum();
        if (scale_val <  h_universe.at(k_dn).at(k)->GetMaximum()) scale_val =  h_universe.at(k_dn).at(k)->GetMaximum();

        h_universe.at(k_up).at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        // FIxed scaling for differential cross section
        if (vars.at(var) != "integrated"){
            // h_universe.at(k_up).at(k)->GetYaxis()->SetRangeUser(0, 0.5e-39);
        }

        TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        if (!single_var){
            leg->AddEntry(h_universe.at(k_up).at(k), Form("%s +1 #sigma", label_pretty.c_str()), "l");
            leg->AddEntry(cv_hist_vec.at(var).at(k),         "CV", "l");
            leg->AddEntry(h_universe.at(k_dn).at(k), Form("%s -1 #sigma", label_pretty.c_str()), "l");
        }
        else {
            leg->AddEntry(h_universe.at(k_up).at(k), Form("%s", label_pretty.c_str()), "l");
            leg->AddEntry(cv_hist_vec.at(var).at(k),         "CV", "l");
        }
        
        leg->Draw();

        bottomPad->cd();
        
        // Up ratio to CV
        TH1D* h_err_up = (TH1D *)h_universe.at(k_up).at(k)->Clone("h_ratio_up");
        h_err_up->Add(cv_hist_vec.at(var).at(k), -1);
        h_err_up->Divide(cv_hist_vec.at(var).at(k));
        h_err_up->SetLineWidth(2);
        h_err_up->SetLineColor(kGreen+2);
        h_err_up->Scale(100);

        h_err_up->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
        
        // Down ratio to CV
        TH1D* h_err_dn = (TH1D *)h_universe.at(k_dn).at(k)->Clone("h_ratio_dn");
        h_err_dn->Add(cv_hist_vec.at(var).at(k), -1);
        h_err_dn->Divide(cv_hist_vec.at(var).at(k));
        h_err_dn->SetLineWidth(2);
        h_err_dn->SetLineColor(kRed+2);
        h_err_dn->Scale(100);

        SetRatioOptions(h_err_up);
        h_err_up->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        h_err_up->GetYaxis()->SetTitle("\% change from CV");
        h_err_up->Draw("hist,same");

        if (!single_var) h_err_dn->Draw("hist,same");

        bottomPad->Update();

        // Draw on the percentage errors manually so they dont overlap
        TLatex* text_up, *text_dn;
        for (int bin = 1; bin < h_err_up->GetNbinsX()+1; bin++){
            
            double bin_up_max = h_err_up->GetBinContent(bin);
            double bin_dn_max = h_err_dn->GetBinContent(bin);

            double shift_up = 0;
            if (bin_up_max >= 0) shift_up = 6;
            else shift_up = -10;

            double shift_dn = 0;
            if (bin_dn_max >= 0) shift_dn = 6;
            else shift_dn = -10;

            if (!single_var){
                if (shift_up < 0 && shift_dn < 0) shift_dn-= 10;
                if (shift_up >= 0 && shift_dn >= 0) shift_up+= 10;
            }
            
            text_up = new TLatex(h_err_up->GetXaxis()->GetBinCenter(bin), bin_up_max+shift_up, Form("%4.1f", h_err_up->GetBinContent(bin)));
            text_dn = new TLatex(h_err_dn->GetXaxis()->GetBinCenter(bin), bin_dn_max+shift_dn, Form("%4.1f", h_err_dn->GetBinContent(bin)));
            text_up->SetTextAlign(21);
            text_up->SetTextColor(kGreen+2);
            text_up->SetTextFont(gStyle->GetTextFont());
            text_up->SetTextSize(0.07);
            text_up->Draw();
            text_dn->SetTextAlign(21);
            text_dn->SetTextColor(kRed+2);
            text_dn->SetTextFont(gStyle->GetTextFont());
            text_dn->SetTextSize(0.07);
            if (!single_var) text_dn->Draw();

        }
        

        FillSysVector(label, var, k, h_err_up, h_err_dn);

        c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(), _util.run_period, label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_err_up;
        delete h_err_dn;
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeDetVar(std::string label, int var, int detvar_index, std::string label_pretty){

    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + vars.at(var));

    std::vector<TH1D*> h_universe;
    std::vector<TH1D*> h_CV;
    
    // Resize to the number cross section types e.g. bkg,eff etc.
    h_universe.resize(xsec_types.size());
    h_CV.resize(xsec_types.size());


    TH1D* h_temp;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < h_universe.size(); k++){
        
        // Get the universe histograms
        _util.GetHist(f_nuexsec, h_temp, Form( "%s/%s/h_run%s_CV_0_%s_%s", label.c_str(), vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), xsec_types.at(k).c_str()));
        h_universe.at(k) = (TH1D*)h_temp->Clone();

        double scale_fact = POT_v.at(k_CV) / POT_v.at(detvar_index);
        h_universe.at(k)->Scale(scale_fact);

        // Get the CV histograms
        _util.GetHist(f_nuexsec, h_CV.at(k), Form( "detvar_CV/%s/h_run%s_CV_0_%s_%s", vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k)->SetLineWidth(2);
        h_universe.at(k)->SetLineColor(kGreen+2);
        h_universe.at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        // Customise
        h_CV.at(k)->SetLineWidth(2);
        h_CV.at(k)->SetLineColor(kBlack);
        h_universe.at(k)->SetLineColor(kGreen+2);
    }

    // Get the covariance matrices
    // CalcMatrices(label, var, h_universe);

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c;
    
    // Now we want to draw them
    for (unsigned int k = 0; k < h_universe.size(); k++){
        c = new TCanvas("c", "c", 500, 500);
        topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        _util.SetTPadOptions(topPad, bottomPad);
        
        h_universe.at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));
        h_universe.at(k)->GetXaxis()->SetTitle("");
        h_universe.at(k)->GetXaxis()->SetLabelSize(0);

        h_universe.at(k)->Draw("hist");
        h_CV.at(k)->Draw("hist,same");

        c->Update();

        double scale_val = h_universe.at(k)->GetMaximum();
        if (scale_val < h_CV.at(k)->GetMaximum()) scale_val = h_CV.at(k)->GetMaximum();

        h_universe.at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        // FIxed scaling for differential cross section
        if (vars.at(var) != "integrated"){
            // h_universe.at(k)->GetYaxis()->SetRangeUser(0, 0.5e-39);
        }

        TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_universe.at(k), label_pretty.c_str(), "l");
        leg->AddEntry(h_CV.at(k), "CV", "l");
        leg->Draw();

        bottomPad->cd();
        
        // Up ratio to CV
        TH1D* h_err = (TH1D *)h_universe.at(k)->Clone("h_ratio");
        h_err->Add(h_CV.at(k), -1);
        h_err->Divide(h_CV.at(k));
        h_err->SetLineWidth(2);
        h_err->SetLineColor(kGreen+2);
        h_err->Scale(100);
        h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
        

        SetRatioOptions(h_err);
        h_err->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        h_err->GetYaxis()->SetTitle("\% change from CV");
        h_err->Draw("hist,same");

        bottomPad->Update();

        // Draw on the percentage errors manually so they dont overlap
        TLatex* text_up, *text_dn;
        for (int bin = 1; bin < h_err->GetNbinsX()+1; bin++){
            
            double bin_up_max = h_err->GetBinContent(bin);

            double shift_up = 0;
            if (bin_up_max >= 0) shift_up = 6;
            else shift_up = -10;
            
            text_up = new TLatex(h_err->GetXaxis()->GetBinCenter(bin), bin_up_max+shift_up, Form("%4.1f", h_err->GetBinContent(bin)));
            text_up->SetTextAlign(21);
            text_up->SetTextColor(kGreen+2);
            text_up->SetTextFont(gStyle->GetTextFont());
            text_up->SetTextSize(0.07);
            text_up->Draw();
        }
        

        // FillSysVector(label, var, k, h_err_up, h_err_dn);

        c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(), _util.run_period, label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_err;
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeMultisim(std::string label, int var, std::string label_pretty, int universes){

    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe; // Universe, <gen/sig/xsec etc>
    std::vector<std::vector<TH1D*>> h_err;

    // Set the X Bins
    std::vector<double> bins;
    if (vars.at(var) == "integrated") bins  = { 0.0, 1.1 };
    else bins = _util.reco_shr_bins;
    // Get the number of bins and the right vector
    int const nbins = bins.size()-1;
    double* edges = &bins[0]; // Cast to an array 

    // Create the Fancier looking 2D histograms
    std::vector<TH2D*> h_universe_2D;
    h_universe_2D.resize(xsec_types.size());
    // Set 2D bins for integrated bins
    if (vars.at(var) == "integrated"){

        h_universe_2D.at(k_xsec_sel)           = new TH2D("h_2D_sel",      "", nbins, edges, 100, 2000, 3000);
        h_universe_2D.at(k_xsec_bkg)           = new TH2D("h_2D_bkg",      "", nbins, edges, 100, 200, 1000);
        h_universe_2D.at(k_xsec_gen)           = new TH2D("h_2D_gen",      "", nbins, edges, 100, 4000, 12000);
        h_universe_2D.at(k_xsec_sig)           = new TH2D("h_2D_sig",      "", nbins, edges, 100, 1000, 3000);
        h_universe_2D.at(k_xsec_eff)           = new TH2D("h_2D_eff",      "", nbins, edges, 100, 0.15, 0.3);
        h_universe_2D.at(k_xsec_ext)           = new TH2D("h_2D_ext",      "", nbins, edges, 10, 0, 10);
        h_universe_2D.at(k_xsec_dirt)          = new TH2D("h_2D_dirt",     "", nbins, edges, 15, 0, 15);
        h_universe_2D.at(k_xsec_data)          = new TH2D("h_2D_data",     "", nbins, edges, 80, 0, 140);
        h_universe_2D.at(k_xsec_mcxsec)        = new TH2D("h_2D_mcxsec",   "", nbins, edges, 100, 0.5e-39, 3.0e-39);
        h_universe_2D.at(k_xsec_dataxsec)      = new TH2D("h_2D_dataxsec", "", nbins, edges, 100, 0.5e-39, 3.0e-39);
        
    }
    // Set 2D bins for other
    else {
        h_universe_2D.at(k_xsec_sel)           = new TH2D("h_2D_sel",      "", nbins, edges, 100, 0, 700);
        h_universe_2D.at(k_xsec_bkg)           = new TH2D("h_2D_bkg",      "", nbins, edges, 100, 0, 250);
        h_universe_2D.at(k_xsec_gen)           = new TH2D("h_2D_gen",      "", nbins, edges, 100, 0, 8000);
        h_universe_2D.at(k_xsec_sig)           = new TH2D("h_2D_sig",      "", nbins, edges, 100, 0, 800);
        h_universe_2D.at(k_xsec_eff)           = new TH2D("h_2D_eff",      "", nbins, edges, 100, 0, 0.5);
        h_universe_2D.at(k_xsec_ext)           = new TH2D("h_2D_ext",      "", nbins, edges, 5, 0, 5);
        h_universe_2D.at(k_xsec_dirt)          = new TH2D("h_2D_dirt",     "", nbins, edges, 10, 0, 10);
        h_universe_2D.at(k_xsec_data)          = new TH2D("h_2D_data",     "", nbins, edges, 100, 0, 50);
        h_universe_2D.at(k_xsec_mcxsec)        = new TH2D("h_2D_mcxsec",   "", nbins, edges, 100, 0, 0.5e-39);
        h_universe_2D.at(k_xsec_dataxsec)      = new TH2D("h_2D_dataxsec", "", nbins, edges, 100, 0, 0.5e-39);
    }
    
    
    // Clone the CV histograms so we can work with them without changing them
    std::vector<TH1D*> cv_hist_vec_clone;
    cv_hist_vec_clone.resize(xsec_types.size());
    
    for (unsigned int h_index = 0; h_index < cv_hist_vec_clone.size(); h_index++){
        cv_hist_vec_clone.at(h_index) = (TH1D*)cv_hist_vec.at(var).at(h_index)->Clone(Form("h_%s_clone", xsec_types.at(h_index).c_str() ));
        // Customise
        cv_hist_vec_clone.at(h_index)->SetLineWidth(2);
        cv_hist_vec_clone.at(h_index)->SetLineColor(kBlack);
    }

    
    // Resize to the number of universes
    h_universe.resize(universes);
    
    // Resize each universe
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        h_universe.at(uni).resize(xsec_types.size());
    }
    
    
    // Get the histograms and customise a bit
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            _util.GetHist(f_nuexsec, h_universe.at(uni).at(k), Form( "%s/%s/h_run%s_%s_%i_%s_%s", label.c_str(), vars.at(var).c_str(), _util.run_period, label.c_str(), uni ,vars.at(var).c_str(), xsec_types.at(k).c_str()));

            // Customise
            h_universe.at(uni).at(k)->SetLineWidth(1);
            h_universe.at(uni).at(k)->SetLineColor(kAzure+5);

            // Loop over the bins and fill the 2D histogram
            for (int bin = 0; bin < h_universe.at(uni).at(k)->GetNbinsX(); bin++){
                h_universe_2D.at(k)->Fill( bins.at(bin), h_universe.at(uni).at(k)->GetBinContent(bin+1));
            }
            
        }
    }

    // Get the covariance, correlation, fractional cov matrices
    CalcMatrices(label, var, h_universe );

    // Create vector of systematic error for each histogram
    std::vector<std::vector<double>> sys_err;
    sys_err.resize(xsec_types.size());
    
    // We now want to get the standard deviation of all universes wrt to the cv
    for (unsigned int i = 0; i < sys_err.size(); i++){
        sys_err.at(i).resize(cv_hist_vec.at(var).at(0)->GetNbinsX(), 0.0);
    }

    // Loop over universes
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over histograms
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

            // Loop over the bins
            for (int bin = 1; bin < cv_hist_vec.at(var).at(0)->GetNbinsX()+1; bin++){
                double uni_x_content = h_universe.at(uni).at(k)->GetBinContent(bin);
                double cv_x_content  = cv_hist_vec.at(var).at(k)->GetBinContent(bin);

                sys_err.at(k).at(bin-1) += ( uni_x_content - cv_x_content) * ( uni_x_content - cv_x_content);
            }
            
        }
    }

    // Loop over the histograms
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        
        // loop over the bins
        for (unsigned int bin = 0; bin < sys_err.at(k).size(); bin ++){
            sys_err.at(k).at(bin) = std::sqrt(sys_err.at(k).at(bin) / h_universe.size());
        }
    }

    // Now we have the systematic error computed we should now set the bin error of the CV clone to the std
    // Loop over universes
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over histograms
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

            // Loop over the bins
            for (int bin = 1; bin < cv_hist_vec.at(var).at(0)->GetNbinsX()+1; bin++){
                cv_hist_vec_clone.at(k)->SetBinError(bin, sys_err.at(k).at(bin-1));
            }
            
        }
    }



    TCanvas *c;
    TPad *topPad;
    TPad *bottomPad;

    // Now we want to draw them
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

        c = new TCanvas("c", "c", 500, 500);
        topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        _util.SetTPadOptions(topPad, bottomPad);

        h_universe.at(0).at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));

        double scale_val = h_universe.at(0).at(k)->GetMaximum();
        
        if (scale_val <  h_universe.at(k_dn).at(k)->GetMaximum()) scale_val =  h_universe.at(k_dn).at(k)->GetMaximum();

        // if (scale_val < h_universe.at(uni).at(k)->GetMaximum()) scale_val = h_universe.at(uni).at(k)->GetMaximum();
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec) h_universe_2D.at(k)->SetTitle(Form("%s", var_labels.at(var).c_str()));
        else if (k == k_xsec_eff) h_universe_2D.at(k)->GetYaxis()->SetTitle("Efficiency");
        else h_universe_2D.at(k)->GetYaxis()->SetTitle("Entries");
        h_universe_2D.at(k)->Draw("colz,same");
        h_universe_2D.at(k)->GetXaxis()->SetTitle("");
        h_universe_2D.at(k)->GetXaxis()->SetLabelSize(0);
        h_universe_2D.at(k)->SetTitle(xsec_types_pretty.at(k).c_str());
        h_universe_2D.at(k)->GetYaxis()->SetTitleSize(0.04);
        h_universe_2D.at(k)->GetYaxis()->SetLabelSize(0.05);



        if (xsec_types.at(k) != "data_xsec") {
            cv_hist_vec_clone.at(k)->SetLineColor(kRed+1);
            cv_hist_vec_clone.at(k)->SetLineWidth(2);
            cv_hist_vec_clone.at(k)->SetLineStyle(7);
            cv_hist_vec_clone.at(k)->SetFillStyle(0);
            cv_hist_vec_clone.at(k)->SetMarkerColor(kRed+1);
            cv_hist_vec_clone.at(k)->Draw("E2,hist,same");
        }
        else { 
            cv_hist_vec_clone.at(k)->SetLineColor(kRed+1);
            cv_hist_vec_clone.at(k)->Draw("E,same");
        }

        h_universe.at(0).at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec){
            // h_universe.at(0).at(k)->GetYaxis()->SetRangeUser(0, 0.5e-39);
        }

        TLegend *leg;
        if (xsec_types.at(k)!= "eff") leg = new TLegend(0.5, 0.65, 0.85, 0.8);
        else leg = new TLegend(0.5, 0.1, 0.85, 0.35);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_universe_2D.at(k), Form("%s", label_pretty.c_str()), "");
        if (xsec_types.at(k) != "data_xsec") leg->AddEntry(cv_hist_vec_clone.at(k),           "CV (Sys Only)", "f");
        else                                 leg->AddEntry(cv_hist_vec_clone.at(k),           "CV (Sys Only)", "le");
        leg->Draw();

        bottomPad->cd();

        // Up percent diff to CV
        TH1D* h_err = (TH1D *)cv_hist_vec_clone.at(k)->Clone("h_err");
        
        // Loop over the bins in the up error, and set the bin content to be the percent difference
        for (int g = 1; g < h_err->GetNbinsX()+1; g++){
            h_err->SetBinContent(g, 100 * h_err->GetBinError(g)/h_err->GetBinContent(g));
        }
        h_err->SetLineWidth(2);
        h_err->SetLineColor(kRed+1);
        h_err->GetYaxis()->SetRangeUser(0, 50);

        h_err->GetXaxis()->SetLabelSize(0.13);
        h_err->GetXaxis()->SetTitleOffset(0.9);
        h_err->GetXaxis()->SetTitleSize(0.13);
        h_err->GetYaxis()->SetLabelSize(0.13);
        h_err->GetYaxis()->SetNdivisions(4, 0, 0, kTRUE);
        h_err->GetYaxis()->SetTitleSize(12);
        h_err->GetYaxis()->SetTitleFont(44);
        h_err->GetYaxis()->CenterTitle();
        h_err->GetYaxis()->SetTitleOffset(1.5);
        h_err->SetTitle(" ");
        h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
        h_err->SetMarkerColor(kBlack);
        h_err->SetLineStyle(1);
        h_err->SetLineColor(kBlack);

        
        h_err->GetYaxis()->SetTitle("\% Uncertainty");
        h_err->SetMarkerSize(4);
        if (k != k_xsec_dirt && k != k_xsec_ext) h_err->Draw("hist, text00");
        gStyle->SetPaintTextFormat("4.1f");

        FillSysVector(label, var, k, h_err, h_err);

        gStyle->SetPalette(56);
        gStyle->SetPalette(kBlueGreenYellow);

        if (xsec_types.at(k) != "data_xsec") _util.Draw_ubooneSim(c, 0.40, 0.915, 0.40, 0.915);
        else _util.Draw_Data_POT(c, Data_POT, 0.50, 0.915, 0.50, 0.915);

        c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_universe_2D.at(k);

    }

        
    
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareCVXSec(){

    std::vector<std::string> error_type = {"stat", "sys", "tot"};

    // Loop over the variables
    for (unsigned int var = 0; var < vars.size(); var++){

        // Loop over the error labels
        for (unsigned int err_lab = 0; err_lab < error_type.size(); err_lab++){

            TH1D* h_dataxsec     = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_temp");
            TH1D* h_dataxsec_tot = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_tot_temp"); // For total uncertainty
            TH1D* h_mcxsec       = (TH1D*) cv_hist_vec.at(var).at(k_xsec_mcxsec)  ->Clone("h_mx_xsec_temp");

            TPad *topPad;
            TPad *bottomPad;

            TCanvas * c = new TCanvas("c", "c", 500, 500);
            topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
            bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
            _util.SetTPadOptions(topPad, bottomPad);

            h_dataxsec->SetLineColor(kBlack);
            h_mcxsec  ->SetLineColor(kRed+2);

            // h_dataxsec->GetYaxis()->SetRangeUser(0, 0.5e-39);
            if (vars.at(var) == "integrated") h_dataxsec->GetYaxis()->SetRangeUser(0.5e-39, 2.5e-39);
            else h_dataxsec->GetYaxis()->SetRangeUser(0.0e-39, 0.5e-39);

            h_dataxsec->GetYaxis()->SetLabelSize(0.04);
            h_dataxsec->GetYaxis()->SetTitleSize(14);
            h_dataxsec->GetYaxis()->SetTitleFont(44);
            h_dataxsec->GetYaxis()->SetTitleOffset(1.5);
            h_dataxsec->GetXaxis()->SetTitle("");
            h_dataxsec->GetXaxis()->SetLabelSize(0);
            h_dataxsec->SetMarkerStyle(20);
            h_dataxsec->SetMarkerSize(0.5);
            h_dataxsec->SetMinimum(0);
            
            // Rewrite the errors for data to sys
            if (error_type.at(err_lab) == "sys"){
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            // Overwrite error to stat + sys
            else {
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec_tot->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin) + v_stat_total.at(var).at(k_xsec_dataxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            
            
            
            h_dataxsec->Draw("E,X0");
            if (error_type.at(err_lab) == "tot"){
                h_dataxsec_tot->Draw("E1,same,X0");
                h_dataxsec->Draw("E1,same,X0");

            } 

            h_mcxsec->Draw("hist,same");
            
            TH1D* h_mcxsec_clone = (TH1D *)h_mcxsec->Clone("h_mc_clone");
            h_mcxsec_clone->SetFillColorAlpha(12, 0.15);
            h_mcxsec_clone->Draw("E2,same");
            

            TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if (error_type.at(err_lab) == "stat")      leg->AddEntry(h_dataxsec, "Data (Stat.)", "le");
            else if (error_type.at(err_lab) == "sys") leg->AddEntry(h_dataxsec, "Data (Sys.)", "le");
            else                                      leg->AddEntry(h_dataxsec, "Data (Stat.+Sys.)", "le");
            leg->AddEntry(h_mcxsec_clone,   "MC (Stat.)", "lf");
            leg->Draw();


            bottomPad->cd();
                
            // The percent difference of mc wrt data
            TH1D* h_err = (TH1D *)h_dataxsec->Clone("h_ratio");
            h_err->Add(h_mcxsec, -1);
            h_err->Divide(h_dataxsec);
            h_err->Scale(100);
            h_err->GetYaxis()->SetTitle("Data - MC / Data [\%]");
            h_err->SetLineWidth(2);
            h_err->SetLineColor(kGreen+2);


            bottomPad->SetGridy(kFALSE);
            
            SetRatioOptions(h_err);
            h_err->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
            h_err->SetLineColor(kBlack);
            h_err->GetYaxis()->SetTitleSize(11);
            h_err->GetYaxis()->SetRangeUser(-100, 100);
            
            h_err->GetYaxis()->SetTitleOffset(2.5);
            h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            if (vars.at(var) == "integrated")  h_err->GetXaxis()->SetLabelSize(0);
            h_err->SetMarkerSize(3.0);
            h_err->Draw("hist,text00");

            // Draw the run period on the plot
            _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

            _util.Draw_Data_POT(c, Data_POT, 0.47, 0.915, 0.47, 0.915);

            
            c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison_%s.pdf", _util.run_period, vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), error_type.at(err_lab).c_str() ));
            delete c;
            delete h_dataxsec;
            delete h_mcxsec;
            delete h_err;

        
        }

    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareCVXSecNoRatio(){

    std::vector<std::string> error_type = {"stat", "sys", "tot"};

    // Loop over the variables
    for (unsigned int var = 0; var < vars.size(); var++){

        // Loop over the error labels
        for (unsigned int err_lab = 0; err_lab < error_type.size(); err_lab++){

            TH1D* h_dataxsec     = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_temp");
            TH1D* h_dataxsec_tot = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_tot_temp"); // For total uncertainty
            TH1D* h_mcxsec       = (TH1D*) cv_hist_vec.at(var).at(k_xsec_mcxsec)  ->Clone("h_mx_xsec_temp");

            TPad *topPad;
            TPad *bottomPad;

            TCanvas * c = new TCanvas("c", "c", 500, 500);

            h_dataxsec->SetLineColor(kBlack);
            h_mcxsec  ->SetLineColor(kRed+2);

            // h_dataxsec->GetYaxis()->SetRangeUser(0, 0.5e-39);
            if (vars.at(var) == "integrated") h_dataxsec->GetYaxis()->SetRangeUser(0.5e-39, 2.5e-39);
            else h_dataxsec->GetYaxis()->SetRangeUser(0.0e-39, 0.5e-39);

            _util.IncreaseLabelSize(h_dataxsec, c);
            if (vars.at(var) == "integrated")h_dataxsec->GetXaxis()->SetLabelSize(0);
            h_dataxsec->GetYaxis()->SetTitleSize(0.04);
            h_dataxsec->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            h_dataxsec->SetMarkerStyle(20);
            h_dataxsec->SetMarkerSize(0.5);
            h_dataxsec->SetMinimum(0);
            
            // Rewrite the errors for data to sys
            if (error_type.at(err_lab) == "sys"){
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1)); // 0.01 is to convert back from a percentage
                }

            }
            // Overwrite error to stat + sys
            else {
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec_tot->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin) + v_stat_total.at(var).at(k_xsec_dataxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            
            
            
            h_dataxsec->Draw("E,X0");
            if (error_type.at(err_lab) == "tot"){
                h_dataxsec_tot->Draw("E1,same,X0");
                h_dataxsec->Draw("E1,same,X0");

            } 

            h_mcxsec->Draw("hist,same");
            
            TH1D* h_mcxsec_clone = (TH1D *)h_mcxsec->Clone("h_mc_clone");
            h_mcxsec_clone->SetFillColorAlpha(12, 0.15);
            h_mcxsec_clone->Draw("E2,same");

            // redraw the data so the data is on top of everything
            h_dataxsec->Draw("E,X0,same");
            if (error_type.at(err_lab) == "tot"){
                h_dataxsec_tot->Draw("E1,same,X0,same");
                h_dataxsec->Draw("E1,same,X0,same");

            }
            

            TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if (error_type.at(err_lab) == "stat")     leg->AddEntry(h_dataxsec, "Data (Stat.)", "le");
            else if (error_type.at(err_lab) == "sys") leg->AddEntry(h_dataxsec, "Data (Sys.)", "le");
            else                                      leg->AddEntry(h_dataxsec, "Data (Stat. + Sys.)", "le");
            leg->AddEntry(h_mcxsec_clone,   "MC (Stat.)", "lf");
            leg->Draw();

            // Draw the run period on the plot
            _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

            _util.Draw_Data_POT(c, Data_POT, 0.52, 0.92, 0.52, 0.92);

            
            c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison_%s_no_ratio.pdf", _util.run_period, vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), error_type.at(err_lab).c_str() ));
            delete c;
            delete h_dataxsec;
            delete h_mcxsec;

        
        }

    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialsePlotCV(){

    // Get the CV histograms. These should stay constant througout the code

    cv_hist_vec.resize(vars.size());
    
    for (unsigned int var = 0; var < vars.size(); var++){
        cv_hist_vec.at(var).resize(xsec_types.size());
    }


    // Loop over the vars
    for (unsigned int var = 0; var < vars.size(); var++){
        
        // Loop over the typrs
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            _util.GetHist(f_nuexsec, cv_hist_vec.at(var).at(k), Form( "CV/%s/h_run%s_CV_0_%s_%s", vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), xsec_types.at(k).c_str()));

            if (cv_hist_vec.at(var).at(k) == NULL) std::cout << "Failed to get the histogram!" << std::endl;

            // Customise
            cv_hist_vec.at(var).at(k)->SetLineWidth(2);
            cv_hist_vec.at(var).at(k)->SetLineColor(kBlack);
        }

        // Create the CV directory and draw the CV
        _util.CreateDirectory("/Systematics/CV/" + vars.at(var) + "/");
    }
    
    TCanvas *c_cv;
    
    // Loop over the vars
    for (unsigned int var = 0; var < vars.size(); var++){
        
        // Loop over the types
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            c_cv = new TCanvas("c", "c", 500, 500);
        
            cv_hist_vec.at(var).at(k)->Draw("hist");
            _util.IncreaseLabelSize(cv_hist_vec.at(var).at(k), c_cv);
        
            TLegend *leg = new TLegend(0.6, 0.8, 0.95, 0.9);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(cv_hist_vec.at(var).at(k),         "CV", "l");
            // leg->Draw();

            // Draw the run period on the plot
            _util.Draw_Run_Period(c_cv, 0.86, 0.92, 0.86, 0.92);

            c_cv->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_%s.pdf", _util.run_period, vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), xsec_types.at(k).c_str()));

            delete c_cv;
            delete leg;
        }
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareVariationXSec(std::string label, int var, std::string label_pretty){

    
    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe;
    
    // Resize to the number of universes
    h_universe.resize(2);
    h_universe.at(k_up).resize(xsec_types.size());
    h_universe.at(k_dn).resize(xsec_types.size());

    // Now get the histograms
    std::string label_up = label + "up";
    std::string label_dn = label + "dn";
    

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        _util.GetHist(f_nuexsec, h_universe.at(k_up).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), _util.run_period, label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        if (k == k_xsec_mcxsec) h_universe.at(k_up).at(k)->SetLineStyle(7);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, h_universe.at(k_dn).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), _util.run_period, label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
        if (k == k_xsec_mcxsec) h_universe.at(k_dn).at(k)->SetLineStyle(7);
    }

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c;
    
    // Now we want to draw them
    c = new TCanvas("c", "c", 500, 500);
    topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    _util.SetTPadOptions(topPad, bottomPad);
    // topPad->SetRightMargin(0.10 );
    // bottomPad->SetRightMargin(0.10 );
    
    h_universe.at(k_up).at(k_xsec_dataxsec)->SetTitle(Form("%s", xsec_types_pretty.at(k_xsec_dataxsec).c_str() ));
    h_universe.at(k_up).at(k_xsec_dataxsec)->GetXaxis()->SetTitle("");
    h_universe.at(k_up).at(k_xsec_dataxsec)->GetXaxis()->SetLabelSize(0);



    h_universe.at(k_up).at(k_xsec_dataxsec)->Draw("hist");
    h_universe.at(k_dn).at(k_xsec_dataxsec)->Draw("hist,same");
    h_universe.at(k_up).at(k_xsec_mcxsec)->Draw("hist,same");
    h_universe.at(k_dn).at(k_xsec_mcxsec)->Draw("hist,same");


    // FIxed scaling for differential cross section
    if (vars.at(var) != "integrated"){
        // h_universe.at(k_up).at(k_xsec_dataxsec)->GetYaxis()->SetRangeUser(0, 0.5e-39);
    }

    TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_universe.at(k_up).at(k_xsec_dataxsec), Form("%s +1 #sigma Data", label_pretty.c_str()), "l");
    leg->AddEntry(h_universe.at(k_up).at(k_xsec_mcxsec),   Form("%s +1 #sigma MC",   label_pretty.c_str()), "l");
    leg->AddEntry(h_universe.at(k_dn).at(k_xsec_dataxsec), Form("%s -1 #sigma Data", label_pretty.c_str()), "l");
    leg->AddEntry(h_universe.at(k_dn).at(k_xsec_mcxsec),   Form("%s -1 #sigma MC",   label_pretty.c_str()), "l");
    leg->Draw();

    bottomPad->cd();
    
    // Up ratio to CV
    TH1D* h_err_up = (TH1D *)h_universe.at(k_up).at(k_xsec_dataxsec)->Clone("h_ratio_up");
    h_err_up->Add(h_universe.at(k_up).at(k_xsec_mcxsec), -1);
    h_err_up->Divide(h_universe.at(k_up).at(k_xsec_dataxsec));
    h_err_up->SetLineWidth(2);
    h_err_up->SetLineColor(kGreen+2);
    h_err_up->Scale(100);
    
    // Down ratio to CV
    TH1D* h_err_dn = (TH1D *)h_universe.at(k_dn).at(k_xsec_dataxsec)->Clone("h_ratio_dn");
    h_err_dn->Add(h_universe.at(k_dn).at(k_xsec_mcxsec), -1);
    h_err_dn->Divide(h_universe.at(k_dn).at(k_xsec_dataxsec));
    h_err_dn->SetLineWidth(2);
    h_err_dn->SetLineColor(kRed+2);
    h_err_dn->Scale(100);

    TH1D* h_err = (TH1D *)cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_ratio");
    h_err->Divide(cv_hist_vec.at(var).at(k_xsec_dataxsec));

    SetRatioOptions(h_err_up);
    h_err_up->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
    h_err_up->GetYaxis()->SetRangeUser(-100, 100);
    h_err_up->GetYaxis()->SetTitle("\% change of Data to MC");
    h_err_up->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
    h_err_up->Draw("hist,same");
    h_err_dn->Draw("hist,same");
    h_err->Draw("hist,same");

    c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_data_mc_comparison.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(), _util.run_period, label.c_str(), vars.at(var).c_str() ));

    delete c;
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CalcMatrices(std::string label, int var, std::vector<std::vector<TH1D*>> h_universe ){

    int n_bins = cv_hist_vec.at(var).at(0)->GetNbinsX();

    TH2D* cov  = new TH2D("h_cov",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the covariance matrix
    TH2D* cor  = new TH2D("h_cor",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the correlation matrix


    // Loop over universes
    std::cout << "Universes: " << h_universe.size() << std::endl;
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over the rows
        for (int row = 1; row < cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetNbinsX()+1; row++){
            
            double uni_row = h_universe.at(uni).at(k_xsec_mcxsec)->GetBinContent(row);
            double cv_row  = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(row);

            // Loop over the columns
            for (int col = 1; col < cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetNbinsX()+1; col++){

                double uni_col = h_universe.at(uni).at(k_xsec_mcxsec)->GetBinContent(col);
                double cv_col  = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(col);
                
                double c = (uni_row - cv_row) * (uni_col - cv_col);

                if (uni != h_universe.size()-1)    cov->SetBinContent(row, col, cov->GetBinContent(row, col) + c ); // Fill with variance 
                else cov->SetBinContent(row, col, (cov->GetBinContent(row, col) + c) / h_universe.size());       // Fill with variance and divide by nuni
            
            }

        }
            
    }


    // ------------ Now calculate the correlation matrix ------------
    double cor_bini;
    // loop over rows
    for (int row = 1; row < cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetNbinsX()+1; row++) {
        
        double cii = cov->GetBinContent(row, row);

        // Loop over columns
        for (int col = 1; col < cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetNbinsX()+1; col++) {
            
            double cjj = cov->GetBinContent(col, col);
            
            double n = sqrt(cii * cjj);

            // Catch Zeros, set to arbitary 1.0
            if (n == 0) cor_bini = 0;
            else cor_bini = cov->GetBinContent(row, col) / n;

            cor->SetBinContent(row, col, cor_bini );
        }
    }

    TH2D *frac_cov = (TH2D*) cov->Clone("h_frac_cov");

    // ------------ Now calculate the fractional covariance matrix ------------
    double setbin;
    // loop over rows
    for (int row = 1; row < cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetNbinsX()+1; row++) {
       
       double cii = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(row);

        // Loop over columns
        for (int col = 1; col < cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetNbinsX()+1; col++) {
            double cjj = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(col);
            double n = cii * cjj;

            // Catch Zeros, set to arbitary 0
            if (n == 0) setbin = 0;
            else setbin = frac_cov->GetBinContent(row, col) / n;

            frac_cov->SetBinContent(row, col, setbin );
        }
    }

    gStyle->SetPalette(kBlueGreenYellow);

    // For now only saving the reco electron cov matrices
    if (var == k_var_reco_el_E){
        // Store the covariance matrices so we can add them and get the total
        // This is a Genie Unisim
        if (label == "RPA"              ||
            label == "CCMEC"            ||
            label == "AxFFCCQE"         ||
            label == "VecFFCCQE"        ||
            label == "DecayAngMEC"      ||
            label == "ThetaDelta2Npi"   ||
            label == "ThetaDelta2NRad"  ||
            label == "RPA_CCQE_Reduced" ||
            label == "NormCCCOH"        ||
            label == "NormNCCOH"){
                h_cov_genie_uni->Add(cov);
                h_cov_tot->Add(cov);
                h_cov_sys->Add(cov);

        }
        else if (label == "Horn_curr" ||
                label == "Horn1_x" ||
                label == "Horn1_y" ||
                label == "Beam_spot" ||
                label == "Horn2_x" ||
                label == "Horn2_y" ||
                label == "Horn_Water" ||
                label == "Beam_shift_x" ||
                label == "Beam_shift_y" ||
                label == "Target_z" ||
                label == "Decay_pipe_Bfield"){
                h_cov_beamline->Add(cov);
                h_cov_tot->Add(cov);
                h_cov_sys->Add(cov);
        }
        else if (label == "weightsGenie"){
            h_cov_genie_multi->Add(cov);
            h_cov_tot->Add(cov);
            h_cov_sys->Add(cov);

        }
        else if (label == "weightsReint"){
            h_cov_reint->Add(cov);
            h_cov_tot->Add(cov);
            h_cov_sys->Add(cov);

        }
        else if (label == "weightsPPFX"){
            h_cov_hp->Add(cov);
            h_cov_tot->Add(cov);
            h_cov_sys->Add(cov);

        }
        else if (label == "Dirt"){
            h_cov_dirt->Add(cov);
            h_cov_tot->Add(cov);
            h_cov_sys->Add(cov);

        }
        else if (label == "POT"){
            h_cov_pot->Add(cov);
            h_cov_tot->Add(cov);
            h_cov_sys->Add(cov);

        }
        else {
            std::cout << "Unknown variation specified: " << label << std::endl;
            return;
        }

    }

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    cov->GetXaxis()->CenterLabels(kTRUE);
    cov->GetYaxis()->CenterLabels(kTRUE);
    cov->GetZaxis()->CenterTitle();
    cov->SetTitle("Covariance Matrix");
    cov->GetZaxis()->SetTitle("Covariance [cm^{4}GeV^{2}]");
    cov->GetZaxis()->SetTitleOffset(1.45);
    cov->Draw("colz");
    _util.IncreaseLabelSize(cov, c);
    _util.Draw_ubooneSim(c, 0.30, 0.915, 0.30, 0.915);
    c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_cov.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str()));

    TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
    cor->GetXaxis()->CenterLabels(kTRUE);
    cor->GetYaxis()->CenterLabels(kTRUE);
    cor->GetZaxis()->CenterTitle();
    cor->GetZaxis()->SetTitle("Correlation");
    cor->GetZaxis()->SetTitleOffset(1.3);
    cor->SetTitle("Correlation Matrix");
    // cor->SetMaximum(1);
    // cor->SetMinimum(0);
    cor->Draw("colz, text00");
    _util.IncreaseLabelSize(cor, c2);
    cor->SetMarkerSize(1.3);
    cor->SetMarkerColor(kRed+1);
    gStyle->SetPaintTextFormat("0.3f");
    _util.Draw_ubooneSim(c2, 0.30, 0.915, 0.30, 0.915);
    c2->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_cor.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str()));

    TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);
    frac_cov->GetXaxis()->CenterLabels(kTRUE);
    frac_cov->GetYaxis()->CenterLabels(kTRUE);
    frac_cov->GetZaxis()->CenterTitle();
    frac_cov->Draw("colz, text00");
    frac_cov->GetZaxis()->SetTitle("Frac. Covariance");
    frac_cov->SetTitle("Fractional Covariance Matrix");
    frac_cov->GetZaxis()->SetTitleOffset(1.52);
    _util.IncreaseLabelSize(frac_cov, c3);
    frac_cov->SetMarkerSize(1.3);
    frac_cov->SetMarkerColor(kRed+1);
    gStyle->SetPaintTextFormat("0.3f");
    _util.Draw_ubooneSim(c3, 0.30, 0.915, 0.30, 0.915);
    c3->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_frac_cov.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str()));

    delete cov;
    delete c;
    delete c2;
    delete cor;
    delete c3;
    delete frac_cov;

}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillSysVector(std::string variation, int var, int type, TH1D *h_up, TH1D *h_dn){

    // This is a Genie Unisim
    if (variation == "RPA"              ||
        variation == "CCMEC"            ||
        variation == "AxFFCCQE"         ||
        variation == "VecFFCCQE"        ||
        variation == "DecayAngMEC"      ||
        variation == "ThetaDelta2Npi"   ||
        variation == "ThetaDelta2NRad"  ||
        variation == "RPA_CCQE_Reduced" ||
        variation == "NormCCCOH"        ||
        variation == "NormNCCOH"){

        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the average error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            av_err += std::abs(h_dn->GetBinContent(bin+1));
            av_err /= std::sqrt(2.0); // sqrt 2 since we are using the covariance matrix formalism (cov matrix with 2 universes) -- error is sqrt diag
            v_genie_uni_total.at(var).at(type).at(bin) += av_err*av_err;
            v_sys_total.at(var).at(type).at(bin)       += av_err*av_err;
            
        }
    }
    // This is a beamline unisim
    else if (
        variation == "Horn_curr" ||
        variation == "Horn1_x" ||
        variation == "Horn1_y" ||
        variation == "Beam_spot" ||
        variation == "Horn2_x" ||
        variation == "Horn2_y" ||
        variation == "Horn_Water" ||
        variation == "Beam_shift_x" ||
        variation == "Beam_shift_y" ||
        variation == "Target_z" ||
        variation == "Decay_pipe_Bfield"){

        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the max error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            av_err += std::abs(h_dn->GetBinContent(bin+1));
            av_err /= std::sqrt(2.0); // sqrt 2 since we are using the covariance matrix formalism (cov matrix with 2 universes) -- error is sqrt diag
            v_beamline_total.at(var).at(type).at(bin) += av_err*av_err;
            v_sys_total.at(var).at(type).at(bin)      += av_err*av_err;
            
        }
    }
    else if (variation == "weightsGenie"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the average error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            // av_err += std::abs(h_dn->GetBinContent(bin+1));
            // av_err /= 2.0;
            v_genie_multi_total.at(var).at(type).at(bin) += av_err*av_err;
            v_sys_total.at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "weightsReint"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the average error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            // av_err += std::abs(h_dn->GetBinContent(bin+1));
            // av_err /= 2.0;
            v_reint_total.at(var).at(type).at(bin) += av_err*av_err;
            v_sys_total.at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "weightsPPFX"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the average error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            // av_err += std::abs(h_dn->GetBinContent(bin+1));
            // av_err /= 2.0;
            v_hp_total.at(var).at(type).at(bin) += av_err*av_err;
            v_sys_total.at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "Dirt"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the average error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            // av_err += std::abs(h_dn->GetBinContent(bin+1));
            // av_err /= 2.0;
            v_dirt_total.at(var).at(type).at(bin) += av_err*av_err;
            v_sys_total.at(var).at(type).at(bin)  += av_err*av_err;
            
        }

    }
    else if (variation == "POT"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the average error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            // av_err += std::abs(h_dn->GetBinContent(bin+1));
            // av_err /= 2.0;
            v_pot_total.at(var).at(type).at(bin) += av_err*av_err;
            v_sys_total.at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else {
        std::cout << "Unknown variation specified: " << variation << std::endl;
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillStatVector(){

    // Loop over the differential variables
    for (unsigned int var = 0; var < vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            
            // Get the uncertainty in each bin
            for (int bin = 0; bin < cv_hist_vec.at(var).at(type)->GetNbinsX(); bin++){
            
                double stat_err = 100 * cv_hist_vec.at(var).at(type)->GetBinError(bin+1) / cv_hist_vec.at(var).at(type)->GetBinContent(bin+1);
                v_stat_total.at(var).at(type).at(bin) += stat_err*stat_err;

            }
        
        }
    }

    // Lets also fill the diagonals of the statistical covariance matrix
    // loop over rows
    for (int row = 1; row < h_cov_stat->GetNbinsX()+1; row++) {
        
        // Loop over columns
        for (int col = 1; col < h_cov_stat->GetNbinsY()+1; col++) {
            
            // Only set the bin content of the diagonals
            double bin_diag = cv_hist_vec.at(k_var_reco_el_E).at(k_xsec_mcxsec)->GetBinContent(row); // We use the MC value to fill the cov matrix, but use the data stat err for now. 

            // 0.01 converts each percentage back to a number. We multiply this by the cv to get the deviate
            if (row == col) h_cov_stat->SetBinContent(row, col, 0.01*0.01*v_stat_total.at(k_var_reco_el_E).at(k_xsec_dataxsec).at(row-1)*bin_diag*bin_diag);  
        }
    }

    // Add the stat to the total
    h_cov_tot->Add(h_cov_stat);


}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillPOTCountingVector(){

    // Loop over the differential variables
    for (unsigned int var = 0; var < vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            
            // Get the uncertainty in each bin
            for (int bin = 0; bin < cv_hist_vec.at(var).at(type)->GetNbinsX(); bin++){
            
                v_pot_total.at(var).at(type).at(bin) += 2.0*2.0;
                v_sys_total.at(var).at(type).at(bin) += 2.0*2.0;

            }
        
        }
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PrintUncertaintySummary(){

    // Loop over the variables
    for (unsigned int var = 0; var < vars.size(); var++ ){
        if (vars.at(var) == "true_el_E") continue; // Skip the true var which doesnt make much sense

        std::cout <<"----------------------------------------------" << std::endl;
        std::cout <<"Differential Variable: " << vars.at(var) <<"\n"<< std::endl;
        
        // Loop over the bins
        for (unsigned int bin = 0; bin < v_genie_uni_total.at(var).at(k_xsec_mcxsec).size(); bin++ ){
            
            std::cout << "Bin: " << bin+1 << " GENIE Unisim:       " <<std::sqrt(v_genie_uni_total.at(var).at(k_xsec_mcxsec)  .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " GENIE Multisim:     " <<std::sqrt(v_genie_multi_total.at(var).at(k_xsec_mcxsec).at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Beamline:           " <<std::sqrt(v_beamline_total.at(var).at(k_xsec_mcxsec)   .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Hadron Prod.:       " <<std::sqrt(v_hp_total.at(var).at(k_xsec_mcxsec)         .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Geant Rein.:        " <<std::sqrt(v_reint_total.at(var).at(k_xsec_mcxsec)      .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Dirt:               " <<std::sqrt(v_dirt_total.at(var).at(k_xsec_mcxsec)       .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " POT Counting:       " <<std::sqrt(v_pot_total.at(var).at(k_xsec_mcxsec)        .at(bin)) << " \%"<< std::endl;
            std::cout << std::endl;
            std::cout << "Bin: " << bin+1 << " Tot Data X-Sec Stat:                 " <<std::sqrt(v_stat_total.at(var).at(k_xsec_dataxsec).at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Stat:                   " <<std::sqrt(v_stat_total.at(var).at(k_xsec_mcxsec).at(bin))   << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Sys:                    " <<std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec)   .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Uncertainty:            " <<std::sqrt(v_stat_total.at(var).at(k_xsec_dataxsec).at(bin) + v_sys_total.at(var).at(k_xsec_mcxsec)   .at(bin)) << " \%"<< std::endl;
            std::cout <<"\n--" << std::endl;       
        }
    }
    std::cout <<"----------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseUncertaintyVectors(){

    // differential variable, type, bin error
    v_genie_uni_total  .resize(vars.size());
    v_genie_multi_total.resize(vars.size());
    v_beamline_total   .resize(vars.size());
    v_hp_total         .resize(vars.size());
    v_reint_total      .resize(vars.size());
    v_sys_total        .resize(vars.size());
    v_stat_total       .resize(vars.size());
    v_dirt_total       .resize(vars.size());
    v_pot_total        .resize(vars.size());

    for (unsigned int var = 0; var < vars.size(); var++ ){
        v_genie_uni_total.at(var)  .resize(xsec_types.size());
        v_genie_multi_total.at(var).resize(xsec_types.size());
        v_beamline_total.at(var)   .resize(xsec_types.size());
        v_hp_total.at(var)         .resize(xsec_types.size());
        v_reint_total.at(var)      .resize(xsec_types.size());
        v_sys_total.at(var)        .resize(xsec_types.size());
        v_stat_total.at(var)       .resize(xsec_types.size());
        v_dirt_total.at(var)       .resize(xsec_types.size());
        v_pot_total .at(var)       .resize(xsec_types.size());
    }

    for (unsigned int var = 0; var < vars.size(); var++ ){
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            v_genie_uni_total.at(var).at(type)  .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_genie_multi_total.at(var).at(type).resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_beamline_total.at(var).at(type)   .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_hp_total.at(var).at(type)         .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_reint_total.at(var).at(type)      .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_sys_total.at(var).at(type)        .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_stat_total.at(var).at(type)       .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_dirt_total.at(var).at(type)       .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_pot_total .at(var).at(type)       .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
        }
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::SaveCovMatrix(TH2D* cov, std::string print_name){

    gStyle->SetPalette(kBlueGreenYellow);

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    cov->GetXaxis()->CenterLabels(kTRUE);
    cov->GetYaxis()->CenterLabels(kTRUE);
    cov->GetZaxis()->CenterTitle();
    cov->GetZaxis()->SetTitle("Covariance [cm^{4}GeV^{2}]");
    cov->GetZaxis()->SetTitleOffset(1.45);
    cov->SetTitle("Covariance Matrix");
    cov->Draw("colz");
    _util.IncreaseLabelSize(cov, c);
    _util.Draw_ubooneSim(c, 0.30, 0.915, 0.30, 0.915);
    c->Print(print_name.c_str());
    delete c;

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotTotUnisim(std::string unisim_type){

    std::vector<std::vector<std::vector<double>>> v_unisim;

    std::vector<std::string> unisim_names;

    // Use the beamline errors
    if (unisim_type == "Beamline")   {
        v_unisim = v_beamline_total;
        unisim_names = {
                    "Horn1_x",
                    "Horn_curr",
                    "Horn1_y",
                    "Beam_spot",
                    "Horn2_x",
                    "Horn2_y",
                    "Horn_Water",
                    "Beam_shift_x",
                    "Beam_shift_y",
                    "Target_z",
                    "Decay_pipe_Bfield"
                };
    }
    else if (unisim_type == "Genie_Unisim") {
        v_unisim = v_genie_uni_total;

        unisim_names = {
                    "RPA",
                    "CCMEC",
                    "AxFFCCQE",
                    "VecFFCCQE",
                    "DecayAngMEC",
                    "ThetaDelta2Npi",
                    "ThetaDelta2NRad",
                    "RPA_CCQE_Reduced",
                    "NormCCCOH",
                    "NormNCCOH"
                };
    }
    else {
        std::cout << "Unknown unisim type specified" << std::endl;
        return;
    }

    // Get all the universes so we can draw them on
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> h_universe; // var -- label -- Up/Dn -- Type
    
    h_universe.resize(vars.size());
    
    for (unsigned int var = 0; var < h_universe.size(); var++){
        h_universe.at(var).resize(unisim_names.size());
    }
    
    // Loop over the vars
    for (unsigned int var = 0; var < h_universe.size(); var++){
        
        // Loop over and resize
        for (unsigned int label = 0; label < h_universe.at(var).size(); label++){
            
            // Resize to the number of universes
            h_universe.at(var).at(label).resize(2);

            // And resize to each type
            h_universe.at(var).at(label).at(k_up).resize(xsec_types.size());
            h_universe.at(var).at(label).at(k_dn).resize(xsec_types.size());
            
        }
    }

    std::vector<std::string> labels_up_v;
    std::vector<std::string> labels_dn_v;

    for (unsigned int var = 0; var < h_universe.size(); var++){

        // Loop over and resize
        for (unsigned int label = 0; label < h_universe.at(var).size(); label++){
        
            // Now get the histograms
            std::string label_up = unisim_names.at(label) + "up";
            std::string label_dn = unisim_names.at(label) + "dn";

            // Set the Unisim up down variation name
            SetLabelName(unisim_names.at(label), label_up, label_dn);

            if (var == 0) labels_up_v.push_back(label_up); // only push back once
            if (var == 0) labels_dn_v.push_back(label_dn); // only push back once


            // Check if its just a single on/off type variation
            // This case we dont want to plot the up/dn, but just once
            bool single_var = false;
            if (label_up == label_dn) single_var = true;

            // Get the histograms and customise a bit
            for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

                TH1D* htemp;
                _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), _util.run_period, label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

                h_universe.at(var).at(label).at(k_up).at(k) = (TH1D*) htemp->Clone(Form("h_clone_%s_%s_up", vars.at(var).c_str(), labels_up_v.at(label).c_str()));

                // Customise
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineWidth(2);
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineStyle(0);
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineColor(kGreen+2);
                
                _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), _util.run_period, label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

                h_universe.at(var).at(label).at(k_dn).at(k) = (TH1D*) htemp->Clone(Form("h_clone_%s_%s_up", vars.at(var).c_str(), labels_dn_v.at(label).c_str()));

                // Customise
                h_universe.at(var).at(label).at(k_dn).at(k)->SetLineWidth(2);
                h_universe.at(var).at(label).at(k_dn).at(k)->SetLineStyle(1);
                h_universe.at(var).at(label).at(k_dn).at(k)->SetLineColor(kRed+2);
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineColor(kGreen+2);
            }
        }
    }


    // Loop over the differential variables
    for (unsigned int var = 0; var < cv_hist_vec.size(); var++ ){
        
        
        // Loop over the types
        for (unsigned int  type = 0; type < cv_hist_vec.at(var).size(); type++ ){

            if (vars.at(var) == "true_el_E" && type != k_xsec_eff) continue; // Skip the true var which doesnt make much sense

            // Get the CV histogram
            TH1D* h_CV_clone = (TH1D*)cv_hist_vec.at(var).at(type)->Clone("h_clone");

            // Loop over the bins and set the error
            for (int bin = 0; bin < h_CV_clone->GetNbinsX(); bin++ ){
                h_CV_clone->SetBinError(bin+1,  0.01*std::sqrt(v_unisim.at(var).at(type).at(bin)) * h_CV_clone->GetBinContent(bin+1));
            }

            // Now plot the damn thing
            TCanvas *c = new TCanvas("c", "c", 500, 500);
            TPad *topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
            TPad *bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
            _util.SetTPadOptions(topPad, bottomPad);
            
            h_CV_clone->SetTitle(Form("%s", unisim_type.c_str() ));

            // if (scale_val < h_universe.at(uni).at(k)->GetMaximum()) scale_val = h_universe.at(uni).at(k)->GetMaximum();
            if (type == k_xsec_mcxsec || type == k_xsec_dataxsec) h_CV_clone->SetTitle(Form("%s", var_labels.at(var).c_str()));
            else if (type == k_xsec_eff) h_CV_clone->GetYaxis()->SetTitle("Efficiency");
            else h_CV_clone->GetYaxis()->SetTitle("Entries");
            h_CV_clone->SetFillColorAlpha(12, 0.15);
            h_CV_clone->Draw("e2, same");
            h_CV_clone->GetXaxis()->SetTitle("");
            h_CV_clone->GetXaxis()->SetLabelSize(0);
            h_CV_clone->SetTitle(xsec_types_pretty.at(type).c_str());
            h_CV_clone->GetYaxis()->SetTitleSize(0.04);
            h_CV_clone->GetYaxis()->SetLabelSize(0.05);

            TLegend *leg;
            
            if (type != k_xsec_eff) leg = new TLegend(0.41, 0.55, 0.91, 0.85);
            else leg = new TLegend(0.4, 0.3, 0.9, 0.6);
            leg->SetNColumns(2);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            leg->AddEntry(h_CV_clone,"CV", "lf");

            // Now draw all the universes
            for (unsigned int label = 0; label < h_universe.at(var).size(); label ++ ){
                SetUnisimColours(unisim_names.at(label), h_universe.at(var).at(label).at(k_up).at(type), h_universe.at(var).at(label).at(k_dn).at(type));
                h_universe.at(var).at(label).at(k_up).at(type)->Draw("hist,same");
                h_universe.at(var).at(label).at(k_dn).at(type)->Draw("hist,same");

                leg->AddEntry(h_universe.at(var).at(label).at(k_up).at(type),Form("%s +1#sigma", unisim_names.at(label).c_str()), "l");
                leg->AddEntry(h_universe.at(var).at(label).at(k_dn).at(type),Form("%s -1#sigma", unisim_names.at(label).c_str()), "l");
            }

            leg->Draw();

            // Draw it again so its on top of everything
            TH1D* h_CV_clone_clone = (TH1D*)h_CV_clone->Clone("h_clone_clone");
            h_CV_clone_clone->SetFillColorAlpha(12, 0.0);
            h_CV_clone_clone->Draw("hist, same");

            bottomPad->cd();

            // Up percent diff to CV
            TH1D* h_err = (TH1D *)h_CV_clone->Clone("h_err");
            
            // Loop over the bins in the up error, and set the bin content to be the percent difference
            for (int g = 1; g < h_err->GetNbinsX()+1; g++){
                h_err->SetBinContent(g, std::sqrt(v_unisim.at(var).at(type).at(g-1)));
            }
            h_err->SetLineWidth(2);
            h_err->SetLineColor(kRed+1);
            h_err->GetYaxis()->SetRangeUser(0, 50);

            h_err->GetXaxis()->SetLabelSize(0.13);
            h_err->GetXaxis()->SetTitleOffset(0.9);
            h_err->GetXaxis()->SetTitleSize(0.13);
            h_err->GetYaxis()->SetLabelSize(0.13);
            h_err->GetYaxis()->SetNdivisions(4, 0, 0, kTRUE);
            h_err->GetYaxis()->SetTitleSize(12);
            h_err->GetYaxis()->SetTitleFont(44);
            h_err->GetYaxis()->CenterTitle();
            h_err->GetYaxis()->SetTitleOffset(1.5);
            h_err->SetTitle(" ");
            h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            h_err->SetMarkerColor(kBlack);
            h_err->SetLineStyle(1);
            h_err->SetLineColor(kBlack);
            h_err->SetFillColorAlpha(12, 0.0);

            
            h_err->GetYaxis()->SetTitle("\% Uncertainty");
            h_err->SetMarkerSize(4);
            if (type != k_xsec_dirt && type != k_xsec_ext) h_err->Draw("hist, text00");
            gStyle->SetPaintTextFormat("4.1f");


            if (xsec_types.at(type) != "data_xsec") _util.Draw_ubooneSim(c, 0.40, 0.915, 0.40, 0.915);
            else _util.Draw_Data_POT(c, Data_POT, 0.50, 0.915, 0.50, 0.915);


            c->Print(Form("plots/run%s/Systematics/%s/run%s_%s_%s_%s.pdf", _util.run_period, unisim_type.c_str(), _util.run_period, unisim_type.c_str(), vars.at(var).c_str(), xsec_types.at(type).c_str()));


            delete c;
            delete h_CV_clone;


        }
    
    
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetUnisimColours(std::string label, TH1D* h_up, TH1D* h_dn){

    if (label == "Horn1_x" || label == "RPA"){
        h_up->SetLineColor(30);
        h_dn->SetLineColor(30);
    }
    else if (label == "Horn_curr" || label == "CCMEC"){
        h_up->SetLineColor(38);
        h_dn->SetLineColor(38);
    }
    else if (label == "Horn1_y" || label == "AxFFCCQE"){
        h_up->SetLineColor(28);
        h_dn->SetLineColor(28);
    }
    else if (label == "Beam_spot" || label == "VecFFCCQE"){
        h_up->SetLineColor(4);
        h_dn->SetLineColor(4);
    }
    else if (label == "Horn2_x" || label == "DecayAngMEC"){
        h_up->SetLineColor(36);
        h_dn->SetLineColor(36);
    }
    else if (label == "Horn2_y" || label == "ThetaDelta2Npi"){
        h_up->SetLineColor(1);
        h_dn->SetLineColor(1);
    }
    else if (label == "Horn_Water" || label == "ThetaDelta2NRad"){
        h_up->SetLineColor(46);
        h_dn->SetLineColor(46);
    }
    else if (label == "Beam_shift_x" || label == "RPA_CCQE_Reduced"){
        h_up->SetLineColor(12);
        h_dn->SetLineColor(12);
    }
    else if (label == "Beam_shift_y" || label == "NormCCCOH"){
        h_up->SetLineColor(kViolet - 7);
        h_dn->SetLineColor(kViolet - 7);
    }
    else if (label == "Target_z" || label == "NormNCCOH"){
        h_up->SetLineColor(kAzure-5);
        h_dn->SetLineColor(kAzure-5);
    }
    else if (label == "Decay_pipe_Bfield"){
        h_up->SetLineColor(kGreen+2);
        h_dn->SetLineColor(kGreen+2);
    }
    else { 
        std::cout << "Unknown label specified: " << label << std::endl;
    }

    h_up->SetLineStyle(1);
    h_dn->SetLineStyle(2);


}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseReweightingModeCut(){
    gStyle->SetOptStat(0);

    // Load in the input file
    // Should we add more protection to this command??
    f_nuexsec = TFile::Open( Form("files/crosssec_run%s.root", _util.run_period ), "READ");

    // Loop over cuts and get the sys uncertainty
    
    for (int cut = 0; cut < _util.k_cuts_MAX ; cut++){
        for (int var = 0; var < _util.k_cut_vars_max; var++){
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "RPA",              2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "CCMEC",            2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "AxFFCCQE",         2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "VecFFCCQE",        2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "DecayAngMEC",      2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "ThetaDelta2Npi",   2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "ThetaDelta2NRad",  2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "RPA_CCQE_Reduced", 2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "NormCCCOH",        2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "NormNCCOH",        2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Horn1_x",          2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Horn_curr",        2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Horn1_y",          2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Beam_spot",        2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Horn2_x",          2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Horn2_y",          2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Horn_Water",       2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Beam_shift_x",     2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Beam_shift_y",     2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Target_z",         2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "Decay_pipe_Bfield",2, "unisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "weightsGenie",600,  "multisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "weightsReint",1000,  "multisim");
            GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, "weightsPPFX", 600, "multisim");
        }
    }
        

}
// -----------------------------------------------------------------------------
void SystematicsHelper::GetCutSysUncertainty(std::string histname, int cut_index, std::string label, int num_uni, std::string var_type){

    // Declare the histogram vector for the cut
    std::vector<TH1D*> h_universe;
    TH1D *h_cv;

    h_universe.resize(num_uni);

    // Now get the histograms
    std::string label_up = label + "up";
    std::string label_dn = label + "dn";

    // Set the Unisim up down variation name
    SetLabelName(label, label_up, label_dn);

    // Now get the histograms

    // If its a unisim then we have up/dn type variation
    if (var_type == "unisim"){
        _util.GetHist(f_nuexsec, h_universe.at(k_up), Form( "%s/Cuts/%s/%s/%s_%s_%s_0", label_up.c_str(), _util.cut_dirs.at(cut_index).c_str(), histname.c_str(), histname.c_str(), label_up.c_str(),_util.cut_dirs.at(cut_index).c_str()));
        _util.GetHist(f_nuexsec, h_universe.at(k_dn), Form( "%s/Cuts/%s/%s/%s_%s_%s_0", label_dn.c_str(), _util.cut_dirs.at(cut_index).c_str(), histname.c_str(), histname.c_str(), label_dn.c_str(),_util.cut_dirs.at(cut_index).c_str()));
    }
    // Multisim
    else {
        // Loop over the universes
        for (int uni = 0; uni < num_uni; uni++){
            _util.GetHist(f_nuexsec, h_universe.at(uni), Form( "%s/Cuts/%s/%s/%s_%s_%s_%i", label.c_str(), _util.cut_dirs.at(cut_index).c_str(), histname.c_str(), histname.c_str(), label.c_str(),_util.cut_dirs.at(cut_index).c_str(), uni));
        }
    }

    // Now get the CV
    _util.GetHist(f_nuexsec, h_cv, Form( "CV/Cuts/%s/%s/%s_CV_%s_0", _util.cut_dirs.at(cut_index).c_str(), histname.c_str(), histname.c_str(),_util.cut_dirs.at(cut_index).c_str()));

    // Now we got the histograms, we loop over an get the uncertainties
    TH1D* h_err = (TH1D*)h_universe.at(k_up)->Clone(); // clone it to get the binning right

    // Clear the bins
    for (int bin = 1; bin < h_universe.at(k_up)->GetNbinsX()+1; bin++){
        h_err->SetBinContent(bin, 0.0);
    }

    // Loop over the universes
    for (unsigned int uni = 0 ; uni < h_universe.size(); uni ++){
        
        // Loop over the bins 
        for (int bin = 1; bin < h_universe.at(k_up)->GetNbinsX()+1; bin++){
            double deviate = h_cv->GetBinContent(bin) - h_universe.at(uni)->GetBinContent(bin); // CV - Uni in bin i
            h_err->SetBinContent(bin, h_err->GetBinContent(bin) + deviate*deviate); // difference squared summed
        }
        
    }

    // Sqrt all bins/N
    for (int bin = 1; bin < h_universe.at(k_up)->GetNbinsX()+1; bin++){
        double err = std::sqrt(h_err->GetBinContent(bin)/num_uni) / h_cv->GetBinContent(bin);
        if (h_cv->GetBinContent(bin) == 0) 
            h_err->SetBinContent(bin, 0.0);
        else 
            h_err->SetBinContent(bin, err);
    }


    // ---- save the histograms into different directories inside the root file
    file_sys_var->cd();

    if(!file_sys_var->GetDirectory(Form("%s/%s", _util.cut_dirs.at(cut_index).c_str(), label.c_str()))) {
        file_sys_var->mkdir(Form("%s/%s", _util.cut_dirs.at(cut_index).c_str(), label.c_str())); // if the directory does not exist, create it
    }

    file_sys_var->cd(Form("%s/%s", _util.cut_dirs.at(cut_index).c_str(), label.c_str())); // open the directory

    h_err->SetDirectory(gDirectory); // set in which dir the hist_ratio.at(k) is going to be written
    h_err->Write(histname.c_str(), TObject::kOverwrite);  // write the histogram to the file

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
