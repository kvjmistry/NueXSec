#include "../include/SystematicsHelper.h"

// -----------------------------------------------------------------------------
void SystematicsHelper::Initialise(Utility _utility){

    std::cout << "Initalising Systematics Helper..." << std::endl;
    _util = _utility;

    // Initialise the Wiener SVD class
    _wSVD.Initialise(_utility);

    // Set the option to scale the bin widths
    if (std::string(_util.scale_bins) == "width")
        scale_bins = true;
    
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
        if (std::string(_util.sysmode) == "default" || std::string(_util.sysmode) == "dedx")  {
            f_vars.at(l) = TFile::Open( Form("files/nuexsec_run%s_%s_merged.root", _util.run_period, var_string.at(l).c_str() ), "READ");
        }
        // Off beam mode
        else if (std::string(_util.sysmode) == "ext") {
            f_vars.at(l) = new TFile( Form("files/nuexsec_ext_run%s_%s.root", _util.run_period, var_string.at(l).c_str() ), "READ");
        }
        // else {
        //     std::cout << "Error I dont know what mode you have configured..." << std::string(_util.sysmode) << std::endl;
        //     exit(1);
        // }
        
        
        // Maybe add something here to pickup non processed variation
    }

    if (std::string(_util.sysmode) == "dedx"){
        MakedEdxPaperPlot();
        return;
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
        // else if (index == k_bnb_diffusion){
        //     h->SetLineColor(kAzure-5);
        // }
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
    
    // Initialise histogram vector to store the uncertainties
    h_cut_err.clear();
    h_cut_err.resize(k_vars_MAX+1); // label -- cut -- variable // plus one to get the total

    for (unsigned int label = 0; label < h_cut_err.size(); label++){
        h_cut_err.at(label).resize(_util.k_cuts_MAX);
    }

    for (unsigned int label = 0; label < h_cut_err.size(); label++){
        
        for (unsigned int cut = 0; cut < h_cut_err.at(label).size(); cut++){
            h_cut_err.at(label).at(cut).resize(_util.k_cut_vars_max);
        }
    }

    for (unsigned int i = 0 ; i < _util.k_cuts_MAX; i++){       
       
        // Default detector systematics mode
        if (mode == "default"){
  
            // Create the directory
            _util.CreateDirectory("/detvar/comparisons/cuts/" + _util.cut_dirs.at(i));

            for(unsigned int j=0; j < _util.vec_hist_name.size(); j++){
                    
                SysVariations(j, Form("plots/run%s/detvar/comparisons/cuts/%s/%s.pdf", _util.run_period, _util.cut_dirs.at(i).c_str(), _util.vec_hist_name.at(j).c_str()), i, _util.vec_axis_label.at(j).c_str(), true);

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

    // Here call a function to write the histograms to file
    if (mode == "default"){
        // Call write function
        SaveCutHistogramsDetVar();
    }


}
// ----------------------------------------------------------------------------
void SystematicsHelper::SysVariations(int hist_index, const char* print_name, int cut, const char* x_axis_name, bool plotdata){

    // ------------------------------------------------------------
    // Some initial configurations and work around fixes
    
    // For some reason the print name keeps changing, this fixes it somewhat
    std::string print_name_str = std::string(print_name);

    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_diff; // The vector of histogram differentes between CV and the vatiation (variation-CV)
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    hist.resize(k_vars_MAX);
    hist_diff.resize(k_vars_MAX);
    hist_ratio.resize(k_vars_MAX);

    // Get the data histogram
    TFile *f_data;
    TH1D *h_data, *h_dirt, *h_ext;

    if (plotdata){
        f_data= TFile::Open("files/nuexsec_run1_merged.root", "READ");
        _util.GetHist(f_data, h_data, Form("Stack/%s/%s/%s_%s_%s", _util.cut_dirs.at(cut).c_str(), "data", _util.vec_hist_name.at(hist_index).c_str(), _util.cut_dirs.at(cut).c_str(), "data"));
        _util.GetHist(f_data, h_dirt, Form("Stack/%s/%s/%s_%s_%s", _util.cut_dirs.at(cut).c_str(), "dirt", _util.vec_hist_name.at(hist_index).c_str(), _util.cut_dirs.at(cut).c_str(), "dirt"));
        _util.GetHist(f_data, h_ext, Form("Stack/%s/%s/%s_%s_%s",  _util.cut_dirs.at(cut).c_str(), "ext",  _util.vec_hist_name.at(hist_index).c_str(), _util.cut_dirs.at(cut).c_str(), "ext"));
        h_data->SetDirectory(0);
        h_dirt->SetDirectory(0);
        h_ext->SetDirectory(0);
        f_data->Close();
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.5);
        h_dirt->Scale(_util.dirt_scale_factor);
        h_ext->Scale(_util.ext_scale_factor);
    }


    TCanvas * c      = new TCanvas("","",500,500);
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    _util.SetTPadOptions(topPad, bottomPad );

    // ------------------------------------------------------------
    // Loop over the variations and get the histograms
    
    for (unsigned int k=0; k < f_vars.size(); k++){
        
        // Loop over the classifications and get the histograms
        for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

            // Only want the MC piece -- may want to add in dirt too? -- will need to separately scale that hitogram though
            if ( i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt ) continue;

            // Get all the MC histograms and add them together
            TH1D *h_temp;

            _util.GetHist(f_vars.at(k), h_temp, Form("Stack/%s/%s/%s_%s_%s", _util.cut_dirs.at(cut).c_str(), _util.classification_dirs.at(i).c_str(),  _util.vec_hist_name.at(hist_index).c_str(), _util.cut_dirs.at(cut).c_str(), _util.classification_dirs.at(i).c_str()));
            
            // First case so clone the histogram
            if (i == 0) hist.at(k) = (TH1D*) h_temp->Clone("h_sum_hist");
            else hist.at(k)->Add(h_temp, 1);
        }
        
    }

    // Legend
    // on top of the topCanvs to avoid overlaping the plot
    TLegend *leg = new TLegend(0.1686747,0.7233083,0.8795181,0.8406015,NULL,"brNDC");
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // ------------------------------------------------------------
    // Now scale the histograms to POT
    // and draw the histograms on the top pad

    // variable used to track the max bin content to later scale the plot to avoid cutting information
    Double_t max_bin = 0.; 

    for (unsigned int y=0; y < hist.size(); y++ ){
        double scale_fact = POT_v.at(k_CV) / POT_v.at(y);
        // std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        if (plotdata){
            hist.at(y)->Scale(_util.mc_scale_factor*3.0754529); // extra factor to scale to main mc pot rather than detvar -- hack
            hist.at(y)->Add(h_dirt, 1.0);
            hist.at(y)->Add(h_ext, 1.0);
        }

        if (y == k_CV){
            // Clone a histogram to plot the CV error as a grey band
            h_error_hist = (TH1D*) hist.at(k_CV)->Clone("h_error_hist");
            h_error_hist->SetFillColorAlpha(12, 0.15);
            h_error_hist->SetLineWidth(2);
            h_error_hist->SetLineColor(kBlack);
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
        if (y == k_CV) {
            leg->AddEntry(h_error_hist, "CV", "lf");
            
            if (hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()) > max_bin)
                max_bin = hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()); // for scale purposes
        
        }
        else {
            hist.at(y)->SetLineColor(var_string_pretty_color.at(y));  // change color of the histogram
            leg->AddEntry(hist.at(y), var_string_pretty.at(y).c_str(), "l"); // add histogram to legend
            
            if (hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()) > max_bin)
                max_bin = hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()); // for scale purposes
        }

    }


    if (plotdata)
         leg->AddEntry(h_data, "Beam-On", "lep"); // add histogram to legend

    // -----------------------------------------------------------------
    // Drawing histograms on top pad

    // setting hist config to the first one that we drawn
    hist.at(0)->GetYaxis()->SetRangeUser(0,max_bin*1.2);
    hist.at(0)->GetYaxis()->SetTitle("Entries");
    hist.at(0)->GetXaxis()->SetLabelSize(0);
    hist.at(0)->GetYaxis()->SetTitleSize(0.05);
    hist.at(0)->GetYaxis()->SetLabelSize(0.05);

    // drawing histograms
    for (unsigned int y=0; y < hist.size(); y++ ) {
    
        if (y == 0)
            hist.at(y)->Draw("hist"); // distinction made so the RangeUser is taken into account, it does not work with "same"
        else
            hist.at(y)->Draw("hist, same");
            //if (y == k_CV) h_error_hist->Draw("E2, same");
    }
    
    // -----------------------------------------------------------------
    // calculate and save the total detector systematics uncertainty 

    // create a temporary histogram that will be used to calculate the total detector sys
    TH1D *h_det_sys_tot;
    
    // loop over variations for a given variable
    for (unsigned int label = 0; label < f_vars.size(); label++){
        
        // ---- save the histograms into different directories inside the root file
        h_cut_err.at(label).at(cut).at(hist_index) = (TH1D*) hist_ratio.at(label)->Clone();
    
        // ---- on the same go, calculate the total detector sys uncertainty

        if( label == 0 ){
            // this is the first histogram, just square it and write it to file
            h_det_sys_tot = (TH1D*) hist_ratio.at(label)->Clone();
            h_det_sys_tot->Multiply(h_det_sys_tot); // square the histogram
        }

        else{
            // this is not the first histogram
            hist_ratio.at(label)->Multiply(hist_ratio.at(label)); // first square the histogram
            h_det_sys_tot->Add(hist_ratio.at(label)); // add it to the existing histogram
        }

    }

    // calculate the square root of the histogram before saving it to the file
    for(int bin = 1; bin <= h_det_sys_tot->GetNbinsX(); bin++){
        h_det_sys_tot->SetBinContent( bin , TMath::Sqrt(h_det_sys_tot->GetBinContent(bin)) );
    }

    // store the total detector sys uncertainty in an additional histogram
    h_cut_err.back().at(cut).at(hist_index) = (TH1D*) h_det_sys_tot->Clone();

    // -----------------------------------------------------------------
    
    // Setting the systematic error
    for(int bin = 1; bin <= h_error_hist->GetNbinsX(); bin++){
        
        // Set the systematic error to be the total det sys error * the bin content
        double sys_err = h_det_sys_tot->GetBinContent(bin) * h_error_hist->GetBinContent(bin);
        h_error_hist->SetBinError(bin, sys_err);
    }

    h_error_hist->Draw("E2, same");
    
    // drawing CV again to make sure it is on top of everything else
    hist.at(k_CV)->Draw("hist, same");

    if (plotdata)
        h_data->Draw("same PE");

    // drawing legend
    leg->Draw();




    // -----------------------------------------------------------------
    // Drawing the total detector systematic uncertainty on the bottom pad

    bottomPad->cd();

    // h_det_sys_tot->GetYaxis()->SetMaxDigits(2);
    h_det_sys_tot->SetLineWidth(2);
    h_det_sys_tot->SetLineColor(1);
    h_det_sys_tot->GetXaxis()->SetLabelSize(15); // 12
    h_det_sys_tot->GetXaxis()->SetLabelFont(43); 
    h_det_sys_tot->GetYaxis()->SetLabelSize(11);
    h_det_sys_tot->GetYaxis()->SetLabelFont(43);
    h_det_sys_tot->GetXaxis()->SetTitleOffset(3.2); // 3
    h_det_sys_tot->GetXaxis()->SetTitleSize(17); // 17
    h_det_sys_tot->GetXaxis()->SetTitleFont(46);
    h_det_sys_tot->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
    h_det_sys_tot->GetYaxis()->SetRangeUser(0, 1.0);
    h_det_sys_tot->GetYaxis()->SetTitle("Tot Frac Unc.");
    h_det_sys_tot->GetYaxis()->SetTitleSize(13); // 13
    h_det_sys_tot->GetYaxis()->SetTitleFont(44);
    h_det_sys_tot->GetYaxis()->SetLabelSize(15); // new
    h_det_sys_tot->GetYaxis()->SetTitleOffset(2);
    h_det_sys_tot->SetTitle(" ");
    h_det_sys_tot->GetXaxis()->SetTitle(x_axis_name);


    h_det_sys_tot->Draw("hist");

    // do you want to print the selection cut stage on your canvas?
    c->cd();
    TLatex *lat = new TLatex(0.15, 0.91, Form("Selection stage: %s", _util.cut_dirs_pretty.at(cut).c_str()));
    lat->SetTextSize(0.03);
    if (!plotdata)lat->Draw();
    c->Modified();

    // Draw the run period on the plot
    if (!plotdata)_util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    if (plotdata)
        _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.45, 0.915, 0.45, 0.915);

    //---------------------------------------------------------------
    // draw final canvas as pdf

    c->Print(print_name_str.c_str()); 

    // close the canvas to avoid warning messages on the terminal
    c->Close();  


    if (plotdata){
        delete h_data;
        delete h_dirt;
        delete h_ext;
    }


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
// ----------------------------------------------------------------------------
void SystematicsHelper::CalcdEdxRMSMean(TH1D* hist, std::string variation){

    // Calculate the Mean and RMS
    double mean_e =0, mean_g = 0; // electron and photon peak
    double rms_e =0, rms_g = 0; // electron and photon peak

    hist->GetXaxis()->SetRange(7, 12);
    mean_e = hist->GetMean();
    rms_e = hist->GetRMS();

    hist->GetXaxis()->SetRange(14, 20);
    mean_g = hist->GetMean();
    rms_g = hist->GetRMS();

    std::cout << variation <<": El Mean: "<< mean_e << " rms: "<< rms_e <<" Photon Mean: "<< mean_g << " rms: "<< rms_g <<  std::endl;

    hist->GetXaxis()->SetRange(1, hist->GetNbinsX());

}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseReweightingMode(){

    gStyle->SetOptStat(0);

    // Load in the input file
    // cross section ratios
    if (std::string(_util.xsec_bin_mode) == "ratio"){
        f_nuexsec = new TFile( Form("files/crosssec_run%s_ratio.root", _util.run_period ), "READ");
    }
    // Normal cross sections
    else {
        f_nuexsec = new TFile( Form("files/crosssec_run%s.root", _util.run_period ), "READ");
    }

    InitialsePlotCV();

    // Now lets initialse the vectors which will store the total uncertainties
    InitialiseUncertaintyVectors();

    // Initialise the covariance matrix vector
    InitialseCovarianceVector();

    // Loop over the cross-section variables
    for (unsigned int var = 0; var <  _util.vars.size(); var++){

        // if (var == 0) continue;
        // if (var != 2) continue;
        // if (var != 0) continue;

        // Comparison plots for data to MC
        CompareVariationXSec("RPA",              var, "RPA" );
        CompareVariationXSec("CCMEC",            var, "CC MEC" );
        CompareVariationXSec("AxFFCCQE",         var, "Ax FF CCQE" );
        CompareVariationXSec("VecFFCCQE",        var, "Vec FF CCQE" );
        CompareVariationXSec("DecayAngMEC",      var, "Decay Ang MEC" );
        CompareVariationXSec("ThetaDelta2Npi",   var, "Theta Delta 2N #pi" );
        CompareVariationXSec("ThetaDelta2NRad",  var, "Theta Delta 2N Rad" );
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
        PlotReweightingModeUnisim("NormCCCOH",        var, "Norm CC COH" );
        PlotReweightingModeUnisim("NormNCCOH",        var, "Norm NC COH" );
        PlotReweightingModeUnisim("xsr_scc_Fv3",      var, "SCC Fv3" );
        PlotReweightingModeUnisim("xsr_scc_Fa3",      var, "SCC Fa3" );


        // Dirt
        PlotReweightingModeUnisim("Dirt",        var, "Dirt" );

        // POT Counting
        PlotReweightingModeUnisim("POT",        var, "POT Count." );

        // PEXT Norm
        PlotReweightingModeUnisim("EXT",        var, "EXT Norm" );

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
        
        // PlotReweightingModeUnisim("Decay_pipe_Bfield",  var, "Decay pipe Bfield" );

        // Detector Variations
        PlotReweightingModeDetVar("LYRayleigh",                         var, k_LYRayleigh,                         var_string_pretty.at(k_LYRayleigh));
        PlotReweightingModeDetVar("LYDown",                             var, k_LYDown,                             var_string_pretty.at(k_LYDown));
        PlotReweightingModeDetVar("SCE",                                var, k_SCE,                                var_string_pretty.at(k_SCE));
        PlotReweightingModeDetVar("Recomb2",                            var, k_Recomb2,                            var_string_pretty.at(k_Recomb2));
        PlotReweightingModeDetVar("WireModX",                           var, k_WireModX,                           var_string_pretty.at(k_WireModX));
        PlotReweightingModeDetVar("WireModYZ",                          var, k_WireModYZ,                          var_string_pretty.at(k_WireModYZ));
        PlotReweightingModeDetVar("WireModThetaXZ",                     var, k_WireModThetaXZ,                     var_string_pretty.at(k_WireModThetaXZ));
        //PlotReweightingModeDetVar("WireModThetaYZ_withSigmaSplines",    var, k_WireModThetaYZ_withSigmaSplines,    var_string_pretty.at(k_WireModThetaYZ_withSigmaSplines));
        
        // PlotReweightingModeDetVar("LYAttenuation",                      var, k_LYAttenuation,                      var_string_pretty.at(k_LYAttenuation));
        // PlotReweightingModeDetVar("WireModThetaYZ_withoutSigmaSplines", var, k_WireModThetaYZ_withoutSigmaSplines, var_string_pretty.at(k_WireModThetaYZ_withoutSigmaSplines));
        // PlotReweightingModeDetVar("WireModdEdX",                        var, k_WireModdEdX,                        var_string_pretty.at(k_WireModdEdX));

        // Plot the multisims
        PlotReweightingModeMultisim("weightsGenie", var,  "GENIE", 500);
        PlotReweightingModeMultisim("weightsReint", var,  "Geant Reinteractions", 1000);
        PlotReweightingModeMultisim("weightsFlux",  var,  "Hadron Production", 500);
        PlotReweightingModeMultisim("MCStats",      var,  "MC Stats", 1000);
        
    }

    // Get the statistical uncertainties
    FillStatVector();

    // Add the smearing covariance matrix into the reco covariance matrix
    AddSmearCovMatrix();

    // Compare the MC and Data X-Section
    // CompareCVXSec();
    CompareCVXSecNoRatio();

    // Save the total covariance matrices
    _util.CreateDirectory("/Systematics/Covariance/" + _util.vars.at(k_var_recoX));
    _util.CreateDirectory("/Systematics/Correlation/" + _util.vars.at(k_var_recoX));
    _util.CreateDirectory("/Systematics/Covariance/" + _util.vars.at(k_var_trueX));
    _util.CreateDirectory("/Systematics/Correlation/" + _util.vars.at(k_var_trueX));
    _util.CreateDirectory("/Systematics/FracCovariance/" + _util.vars.at(k_var_trueX));
    _util.CreateDirectory("/Systematics/FracCovariance/" + _util.vars.at(k_var_recoX));
    for (unsigned int cov = 0; cov < h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).size(); cov++){
        
        SaveCovMatrix(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(cov),                                                        Form("plots/run%s/Systematics/Covariance/%s/run%s_%s_%s_%s_cov.pdf",          _util.run_period, _util.vars.at(k_var_recoX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_mcxsec).c_str(),   _util.vars.at(k_var_recoX).c_str()));
        SaveCovMatrix(h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(cov),                                                      Form("plots/run%s/Systematics/Covariance/%s/run%s_%s_%s_%s_cov.pdf",          _util.run_period, _util.vars.at(k_var_recoX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_dataxsec).c_str(), _util.vars.at(k_var_recoX).c_str()));
        SaveCorMatrix(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(cov),       cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec),   Form("plots/run%s/Systematics/Correlation/%s/run%s_%s_%s_%s_cor.pdf",         _util.run_period, _util.vars.at(k_var_recoX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_mcxsec).c_str(),   _util.vars.at(k_var_recoX).c_str()));
        SaveCorMatrix(h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(cov),     cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec), Form("plots/run%s/Systematics/Correlation/%s/run%s_%s_%s_%s_cor.pdf",         _util.run_period, _util.vars.at(k_var_recoX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_dataxsec).c_str(), _util.vars.at(k_var_recoX).c_str()));
        SaveFracCovMatrix(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(cov),   cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec),   Form("plots/run%s/Systematics/FracCovariance/%s/run%s_%s_%s_%s_frac_cov.pdf", _util.run_period, _util.vars.at(k_var_recoX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_mcxsec).c_str(),   _util.vars.at(k_var_recoX).c_str()));
        SaveFracCovMatrix(h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(cov), cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec), Form("plots/run%s/Systematics/FracCovariance/%s/run%s_%s_%s_%s_frac_cov.pdf", _util.run_period, _util.vars.at(k_var_recoX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_dataxsec).c_str(), _util.vars.at(k_var_recoX).c_str()));
    
        SaveCovMatrix(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(cov),                                                          Form("plots/run%s/Systematics/Covariance/%s/run%s_%s_%s_%s_cov_shape.pdf",          _util.run_period, _util.vars.at(k_var_trueX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_mcxsec_shape).c_str(),   _util.vars.at(k_var_trueX).c_str()));
        SaveCorMatrix(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(cov),     cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape), Form("plots/run%s/Systematics/Correlation/%s/run%s_%s_%s_%s_cor_shape.pdf",         _util.run_period, _util.vars.at(k_var_trueX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_mcxsec_shape).c_str(), _util.vars.at(k_var_trueX).c_str()));
        SaveFracCovMatrix(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(cov),   cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape),   Form("plots/run%s/Systematics/FracCovariance/%s/run%s_%s_%s_%s_frac_cov_shape.pdf", _util.run_period, _util.vars.at(k_var_trueX).c_str(), _util.run_period, systematic_names.at(cov).c_str(), xsec_types.at(k_xsec_mcxsec_shape).c_str(),   _util.vars.at(k_var_trueX).c_str()));
    }


    // Create the directories
    _util.CreateDirectory("/Systematics/Beamline");
    _util.CreateDirectory("/Systematics/Genie_Unisim");
    _util.CreateDirectory("/Systematics/DetVar");
    
    // Plot the total beamline sys uncertainty
    PlotTotUnisim("Beamline");
    PlotTotUnisim("Genie_Unisim");
    PlotTotUnisim("DetVar");

    // Print a summary of the results
    PrintUncertaintySummary();

    // Plot the systematic uncertainty
    MakeTotUncertaintyPlot(false);
    MakeTotUncertaintyPlot(true);

    // Calculate a chi-squared using the total covariance matrix
    double chi, pval;
    int ndof;
    _util.CalcChiSquared(cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec), cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec), h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_tot), chi, ndof, pval);
    _util.CalcChiSquared(cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec), cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec), h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(k_err_tot), chi, ndof, pval);

    // Export the results to file
    ExportResult(f_nuexsec);
    ExportTotalCrossSectionResult();


}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetLabelName(std::string label, std::string &label_up, std::string &label_dn){

    if  (label == "Horn_curr"          ){
        label_up = "Horn_p2kA";
        label_dn = "Horn_m2kA";
    }
    else if  (label == "Horn1_x"       ){
        label_up = label + "_p3mm";
        label_dn = label + "_m3mm";
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
        label_dn = label + "_m3mm";
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
    else if  (label == "Horn1_refined_descr"                || label == "Decay_pipe_Bfield" || label == "Old_Horn_Geometry" ||
              label == "LYRayleigh"                         || label == "LYDown"            || label == "LYAttenuation"     || label == "SCE" || label == "Recomb2"       || 
              label == "WireModX"                           || label == "WireModYZ"         || label == "WireModThetaXZ"    || label == "WireModThetaYZ_withSigmaSplines" || 
              label == "WireModThetaYZ_withoutSigmaSplines" || label == "WireModdEdX"       || label == "EXT"){
        label_up = label;
        label_dn = label;
    }
    else if  (label == "RPA"           ){
        label_up = label + "up";
        label_dn = label + "dn";
    }
    else if  (label == "CCMEC"         ){
        label_up = label + "dn";
        label_dn = label + "dn";
    }
    else if  (label == "AxFFCCQE"      ){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "VecFFCCQE"     ){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "DecayAngMEC"   ){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "ThetaDelta2Npi"){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "ThetaDelta2NRad"){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "NormCCCOH"      ){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "NormNCCOH"      ){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "xsr_scc_Fv3"      ){
        label_up = label + "up";
        label_dn = label + "up";
    }
    else if  (label == "xsr_scc_Fa3"      ){
        label_up = label + "up";
        label_dn = label + "up";
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
    _util.CreateDirectory("/Systematics/" + label + "/" + _util.vars.at(var));

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

    TH1D* htemp;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), _util.vars.at(var).c_str(), _util.run_period, label_up.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

        h_universe.at(k_up).at(k) = (TH1D*)htemp->Clone();

        if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
            h_universe.at(k_up).at(k)->Scale(1.0, "width");

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), _util.vars.at(var).c_str(), _util.run_period, label_dn.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

        h_universe.at(k_dn).at(k) = (TH1D*)htemp->Clone();

        if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
            h_universe.at(k_dn).at(k)->Scale(1.0, "width");

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
    }

    std::vector<std::vector<TH1D*>> h_universe_uni;
    
    // Get the covariance matrices
    
    // Single var so only 1 universe which we store in "up"
    if (single_var){
        h_universe_uni.push_back(h_universe.at(k_up));
        for (unsigned int type = 0; type < xsec_types.size(); type++) { 
            CalcMatrices(label, var, h_universe_uni, type, cv_hist_vec.at(var).at(type) );
        }
    }
    // Double var, so use up and down
    else {
        for (unsigned int type = 0; type < xsec_types.size(); type++) { 
            CalcMatrices(label, var, h_universe, type, cv_hist_vec.at(var).at(type) );
        }
    }
    
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

        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_universe.at(k_up).at(k)->SetTitle(_util.var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_mcxsec_smear || k == k_xsec_mcxsec_shape)
            h_universe.at(k_up).at(k)->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
        else if (k == k_xsec_eff)
            h_universe.at(k_up).at(k)->SetTitle(_util.var_labels_eff.at(var).c_str());
        else
            h_universe.at(k_up).at(k)->SetTitle(_util.var_labels_events.at(var).c_str());
        
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
        if (_util.vars.at(var) != "integrated"){
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
        
        // Down ratio to CV
        TH1D* h_err_dn = (TH1D *)h_universe.at(k_dn).at(k)->Clone("h_ratio_dn");
        h_err_dn->Add(cv_hist_vec.at(var).at(k), -1);

        // Store the uncertainty
        FillSysVector(label, var, k, h_err_up, h_err_dn, cv_hist_vec.at(var).at(k));

        h_err_up->Divide(cv_hist_vec.at(var).at(k));
        h_err_up->SetLineWidth(2);
        h_err_up->SetLineColor(kGreen+2);
        h_err_up->Scale(100);

        h_err_dn->Divide(cv_hist_vec.at(var).at(k));
        h_err_dn->SetLineWidth(2);
        h_err_dn->SetLineColor(kRed+2);
        h_err_dn->Scale(100);

        

        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_err_up->SetTitle(_util.var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_mcxsec_smear || k == k_xsec_mcxsec_shape)
            h_err_up->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
        else if (k == k_xsec_eff)
            h_err_up->SetTitle(_util.var_labels_eff.at(var).c_str());
        else
            h_err_up->SetTitle(_util.var_labels_events.at(var).c_str());        
        
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
        
        // Choose what histograms to save
        if ( (var == k_var_recoX && (k == k_xsec_mcxsec || k == k_xsec_dataxsec || k == k_xsec_bkg || k == k_xsec_sig)) || 
             (var == k_var_trueX && (k == k_xsec_eff || k == k_xsec_mcxsec_smear || k == k_xsec_mcxsec_shape)) ||
             (var == k_var_integrated && (k == k_xsec_mcxsec || k == k_xsec_dataxsec || k == k_xsec_bkg || k == k_xsec_sig || k == k_xsec_eff)) 
            ){
            c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), _util.vars.at(var).c_str(), _util.run_period, label.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));
        }

        delete c;
        delete leg;
        delete h_err_up;
        delete h_err_dn;
        delete h_universe.at(k_up).at(k);
        delete h_universe.at(k_dn).at(k);
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeDetVar(std::string label, int var, int detvar_index, std::string label_pretty){

    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + _util.vars.at(var));

    std::vector<TH1D*> h_universe;
    std::vector<TH1D*> h_CV;
    
    // Resize to the number cross section types e.g. bkg,eff etc.
    h_universe.resize(xsec_types.size());
    h_CV.resize(xsec_types.size());


    TH1D* h_temp;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < h_universe.size(); k++){
        
        // Get the universe histograms
        _util.GetHist(f_nuexsec, h_temp, Form( "%s/%s/h_run%s_CV_0_%s_%s", label.c_str(), _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));
        h_universe.at(k) = (TH1D*)h_temp->Clone();

        if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
            h_universe.at(k)->Scale(1.0, "width");

        double scale_fact = POT_v.at(k_CV) / POT_v.at(detvar_index);
        
        // Scale the histograms, but only in the case of certain variables
        if (xsec_types.at(k) == "sel" || xsec_types.at(k) == "gen" || xsec_types.at(k) == "bkg" || xsec_types.at(k) == "sig") 
            h_universe.at(k)->Scale(scale_fact);

        // Get the CV histograms
        _util.GetHist(f_nuexsec, h_temp, Form( "detvar_CV/%s/h_run%s_CV_0_%s_%s", _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));
        h_CV.at(k) = (TH1D*)h_temp->Clone();

        if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
            h_CV.at(k)->Scale(1.0, "width");

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

    // In the case of a true variation, do the CV background subtraction 
    if (var == k_var_trueX){
        
        TH1D* h_bkg_CV, *h_ext_CV, *h_dirt_CV;
        
        // Load in the backgrounds from file
        TH1D* h_temp = (TH1D*)f_nuexsec->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_bkg", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
        h_bkg_CV = (TH1D*)h_temp->Clone();

        h_temp = (TH1D*)f_nuexsec->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_ext", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
        h_ext_CV = (TH1D*) h_temp->Clone();
        h_ext_CV->Scale(_util.ext_scale_factor /  (_util.config_v.at(_util.k_Run1_Data_POT)  / POT_v.at(k_CV)) );
        h_bkg_CV->Add(h_ext_CV);
        
        h_temp = (TH1D*)f_nuexsec->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_dirt", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
        h_dirt_CV = (TH1D*)h_temp->Clone();
        h_dirt_CV->Scale(_util.dirt_scale_factor / (_util.config_v.at(_util.k_Run1_Data_POT)  / POT_v.at(k_CV)) );
        h_bkg_CV->Add(h_dirt_CV);

        // HARDCODED -- NEED TO CONFIGURE A GOLBAL CALC OF FLUX AND NUM TARG
        // Scale by the flux and number of targets for the detvar CV
        h_bkg_CV->Scale(1.0/(1.40254e+10 * 4.31247e+31));
        h_bkg_CV->Scale(1.0e39);
        
        if (scale_bins)
            h_bkg_CV->Scale(1.0, "width");

        h_universe.at(k_xsec_mcxsec_shape)->Add(h_bkg_CV, -1);
        delete h_bkg_CV;
    }

    // Get the covariance matrices
    std::vector<std::vector<TH1D*>> h_universe_uni;
    h_universe_uni.push_back(h_universe);
    for (unsigned int type = 0; type < xsec_types.size(); type++) { 
        CalcMatrices(label, var, h_universe_uni, type, h_CV.at(type));
    }

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
        if (_util.vars.at(var) != "integrated"){
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

        // Store the sys error
        FillSysVector(label, var, k, h_err, h_err, h_CV.at(k));

        h_err->Divide(h_CV.at(k));
        h_err->SetLineWidth(2);
        h_err->SetLineColor(kGreen+2);
        h_err->Scale(100);
        
        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_err->SetTitle(_util.var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_mcxsec_smear || k == k_xsec_mcxsec_shape)
            h_err->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
        else if (k == k_xsec_eff)
            h_err->SetTitle(_util.var_labels_eff.at(var).c_str());
        else
            h_err->SetTitle(_util.var_labels_events.at(var).c_str());
        
        // h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
        

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
        
        if ( (var == k_var_recoX && (k == k_xsec_mcxsec || k == k_xsec_dataxsec || k == k_xsec_bkg || k == k_xsec_sig || k == k_xsec_sel)) || 
             (var == k_var_trueX && (k == k_xsec_eff || k == k_xsec_mcxsec_smear || k == k_xsec_mcxsec_shape)) ||
             (var == k_var_integrated && (k == k_xsec_mcxsec || k == k_xsec_dataxsec || k == k_xsec_bkg || k == k_xsec_sig || k == k_xsec_eff)) 
            ){
            c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), _util.vars.at(var).c_str(), _util.run_period, label.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));
        }

        delete c;
        delete leg;
        delete h_err;
        delete h_universe.at(k);
        delete h_CV.at(k);
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeMultisim(std::string label, int var, std::string label_pretty, int universes){

    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + _util.vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe; // Universe, <gen/sig/xsec etc>
    std::vector<std::vector<TH1D*>> h_err;

    // Set the X Bins
    std::vector<double> bins;
    if
        (_util.vars.at(var) == "integrated") bins  = { 0.0, 1.1 };
    else {
        // Electron/Shower Energy
        if (std::string(_util.xsec_var) =="elec_E"){
            bins = _util.reco_shr_bins;
        }
        // Electron/Shower beta
        else if (std::string(_util.xsec_var) =="elec_ang"){
            bins = _util.reco_shr_bins_ang;
        }
        // Electron/Shower cos beta
        else if (std::string(_util.xsec_var) =="elec_cang"){
            bins = _util.reco_shr_bins_cang;
        }
        else {
            std::cout << "Unsupported parameter...exiting!" << std::endl;
            return;
        }
        
    }
    
    // Get the number of bins and the right vector
    int const nbins = bins.size()-1;
    double* edges = &bins[0]; // Cast to an array 

    // Create the Fancier looking 2D histograms
    std::vector<TH2D*> h_universe_2D;
    h_universe_2D.resize(xsec_types.size());
    // Set 2D bins for integrated bins
    if (_util.vars.at(var) == "integrated"){

        std::vector<double> int_bins_low;
        std::vector<double> int_bins_high;
        
        if (std::string(_util.xsec_smear_mode) == "mcc8" ){
            //               sel - bkg -   gen - gen_smear - sig  - eff - ext - dirt - data - mcxsec- mcxsec_smear - dataxsec
            int_bins_low  = {2000, 200,   4000,       1000,  1000, 0.15,   0,     0,      0,     0.5,          0.5,0.5,     0.5}; 
            int_bins_high = {3000, 1000, 12000,       3000,  3000, 0.3,    10,    15,   140,     3.0,          3.0,3.0,     3.0}; 
        }
        // Event Rate Binning
        else {
             //               sel - bkg -   gen - gen_smear - sig  - eff - ext - dirt - data - mcxsec - mcxsec_smear - dataxsec
            int_bins_low  = {3000, 200,   4000,       1000,  1000, 0.1,    0,     0,      0,     0.5,            0.5,0.5,    0.5}; 
            int_bins_high = {4000, 1000, 16000,       3000,  3000, 0.3,    10,    15,   140,     3.0,            3.0,3.0,    3.0}; 

        }

        h_universe_2D.at(k_xsec_sel)           = new TH2D("h_2D_sel",          "", nbins, edges, 100, int_bins_low.at(k_xsec_sel),          int_bins_high.at(k_xsec_sel));
        h_universe_2D.at(k_xsec_bkg)           = new TH2D("h_2D_bkg",          "", nbins, edges, 100, int_bins_low.at(k_xsec_bkg),          int_bins_high.at(k_xsec_bkg));
        h_universe_2D.at(k_xsec_gen)           = new TH2D("h_2D_gen",          "", nbins, edges, 100, int_bins_low.at(k_xsec_gen),          int_bins_high.at(k_xsec_gen));
        h_universe_2D.at(k_xsec_gen_smear)     = new TH2D("h_2D_gen_smear",    "", nbins, edges, 100, int_bins_low.at(k_xsec_gen_smear),    int_bins_high.at(k_xsec_gen_smear));
        h_universe_2D.at(k_xsec_sig)           = new TH2D("h_2D_sig",          "", nbins, edges, 100, int_bins_low.at(k_xsec_sig),          int_bins_high.at(k_xsec_sig));
        h_universe_2D.at(k_xsec_eff)           = new TH2D("h_2D_eff",          "", nbins, edges, 100, int_bins_low.at(k_xsec_eff),          int_bins_high.at(k_xsec_eff));
        h_universe_2D.at(k_xsec_ext)           = new TH2D("h_2D_ext",          "", nbins, edges, 10,  int_bins_low.at(k_xsec_ext),          int_bins_high.at(k_xsec_ext));
        h_universe_2D.at(k_xsec_dirt)          = new TH2D("h_2D_dirt",         "", nbins, edges, 15,  int_bins_low.at(k_xsec_dirt),         int_bins_high.at(k_xsec_dirt));
        h_universe_2D.at(k_xsec_data)          = new TH2D("h_2D_data",         "", nbins, edges, 80,  int_bins_low.at(k_xsec_data),         int_bins_high.at(k_xsec_data));
        h_universe_2D.at(k_xsec_mcxsec)        = new TH2D("h_2D_mcxsec",       "", nbins, edges, 100, int_bins_low.at(k_xsec_mcxsec),       int_bins_high.at(k_xsec_mcxsec));
        h_universe_2D.at(k_xsec_mcxsec_smear)  = new TH2D("h_2D_mcxsec_smear", "", nbins, edges, 100, int_bins_low.at(k_xsec_mcxsec_smear), int_bins_high.at(k_xsec_mcxsec_smear));
        h_universe_2D.at(k_xsec_mcxsec_shape)  = new TH2D("h_2D_mcxsec_shape", "", nbins, edges, 100, int_bins_low.at(k_xsec_mcxsec_shape), int_bins_high.at(k_xsec_mcxsec_shape));
        h_universe_2D.at(k_xsec_dataxsec)      = new TH2D("h_2D_dataxsec",     "", nbins, edges, 100, int_bins_low.at(k_xsec_dataxsec),     int_bins_high.at(k_xsec_dataxsec));
        
    }
    // Set 2D bins for other
    else {

        std::vector<double> diff_bins_low;
        std::vector<double> diff_bins_high;
        
        if (std::string(_util.xsec_smear_mode) == "mcc8" ){
            //               sel - bkg -   gen - gen_smear - sig  - eff - ext - dirt - data - mcxsec - mcxsec_smear/shape - dataxsec
            diff_bins_low  = { 0,    0,      0,          0,    0,     0,    0,     0,      0,       0,            0,0,     0}; 
            diff_bins_high = {2500, 1200, 16000,      2700, 2700,     5,   30,    40,    250, _util.xsec_scale, _util.xsec_scale, _util.xsec_scale,  _util.xsec_scale}; 
        }
        // Event Rate Binning
        else {
            //               sel - bkg -   gen - gen_smear - sig  - eff - ext - dirt - data - mcxsec - mcxsec_smear/shape - dataxsec
            diff_bins_low  = { 0,    0,      0,          0,    0,     0,    0,     0,      0,      0,             0,0,        0}; 
            diff_bins_high = {4000, 1200, 16000,      2700, 2700,     0.5, 50,    50,    250, _util.xsec_scale, _util.xsec_scale, _util.xsec_scale, _util.xsec_scale}; 

        }

        h_universe_2D.at(k_xsec_sel)           = new TH2D("h_2D_sel",          "", nbins, edges, 100, diff_bins_low.at(k_xsec_sel),           diff_bins_high.at(k_xsec_sel));
        h_universe_2D.at(k_xsec_bkg)           = new TH2D("h_2D_bkg",          "", nbins, edges, 100, diff_bins_low.at(k_xsec_bkg),           diff_bins_high.at(k_xsec_bkg));
        h_universe_2D.at(k_xsec_gen)           = new TH2D("h_2D_gen",          "", nbins, edges, 100, diff_bins_low.at(k_xsec_gen),           diff_bins_high.at(k_xsec_gen));
        h_universe_2D.at(k_xsec_gen_smear)     = new TH2D("h_2D_gen_smear",    "", nbins, edges, 100, diff_bins_low.at(k_xsec_gen_smear),     diff_bins_high.at(k_xsec_gen_smear));
        h_universe_2D.at(k_xsec_sig)           = new TH2D("h_2D_sig",          "", nbins, edges, 100, diff_bins_low.at(k_xsec_sig),           diff_bins_high.at(k_xsec_sig));
        h_universe_2D.at(k_xsec_eff)           = new TH2D("h_2D_eff",          "", nbins, edges, 100, diff_bins_low.at(k_xsec_eff),           diff_bins_high.at(k_xsec_eff));
        h_universe_2D.at(k_xsec_ext)           = new TH2D("h_2D_ext",          "", nbins, edges, 30,  diff_bins_low.at(k_xsec_ext),           diff_bins_high.at(k_xsec_ext));
        h_universe_2D.at(k_xsec_dirt)          = new TH2D("h_2D_dirt",         "", nbins, edges, 40,  diff_bins_low.at(k_xsec_dirt),          diff_bins_high.at(k_xsec_dirt));
        h_universe_2D.at(k_xsec_data)          = new TH2D("h_2D_data",         "", nbins, edges, 100, diff_bins_low.at(k_xsec_data),          diff_bins_high.at(k_xsec_data));
        h_universe_2D.at(k_xsec_mcxsec)        = new TH2D("h_2D_mcxsec",       "", nbins, edges, 100, diff_bins_low.at(k_xsec_mcxsec),        diff_bins_high.at(k_xsec_mcxsec));
        h_universe_2D.at(k_xsec_mcxsec_smear)  = new TH2D("h_2D_mcxsec_smear", "", nbins, edges, 100, diff_bins_low.at(k_xsec_mcxsec_smear),  diff_bins_high.at(k_xsec_mcxsec_smear));
        h_universe_2D.at(k_xsec_mcxsec_shape)  = new TH2D("h_2D_mcxsec_shape", "", nbins, edges, 100, diff_bins_low.at(k_xsec_mcxsec_shape),  diff_bins_high.at(k_xsec_mcxsec_shape));
        h_universe_2D.at(k_xsec_dataxsec)      = new TH2D("h_2D_dataxsec",     "", nbins, edges, 100, diff_bins_low.at(k_xsec_dataxsec),      diff_bins_high.at(k_xsec_dataxsec));
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
    for (int uni = 0; uni < universes; uni++){
        h_universe.at(uni).resize(xsec_types.size());
    }
    
    
    // Get the histograms and customise a bit
    for (int uni = 0; uni < universes; uni++){
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            _util.GetHist(f_nuexsec, h_universe.at(uni).at(k), Form( "%s/%s/h_run%s_%s_%i_%s_%s", label.c_str(), _util.vars.at(var).c_str(), _util.run_period, label.c_str(), uni ,_util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

            if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
                h_universe.at(uni).at(k)->Scale(1.0, "width");

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
    std::cout << "Calculating Covariance Matrix" << std::endl;
    for (unsigned int type = 0; type < xsec_types.size(); type++) {  
        CalcMatrices(label, var, h_universe, type, cv_hist_vec.at(var).at(type));
    }
    std::cout << "Finished Calculating Covariance Matrix" << std::endl;

    // Create vector of systematic error for each histogram
    std::vector<std::vector<double>> sys_err;
    sys_err.resize(xsec_types.size());
    
    // We now want to get the standard deviation of all universes wrt to the cv
    for (unsigned int i = 0; i < sys_err.size(); i++){
        sys_err.at(i).resize(cv_hist_vec.at(var).at(0)->GetNbinsX(), 0.0);
    }

    // Loop over universes
    for (int uni = 0; uni < universes; uni++){
        
        // Loop over histograms
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

            // Loop over the bins
            for (int bin = 1; bin < cv_hist_vec.at(var).at(0)->GetNbinsX()+1; bin++){
                double uni_x_content = h_universe.at(uni).at(k)->GetBinContent(bin);
                double cv_x_content  = cv_hist_vec.at(var).at(k)->GetBinContent(bin);

                sys_err.at(k).at(bin-1) += ( uni_x_content - cv_x_content) * ( uni_x_content - cv_x_content);
            }
            
            // Done with the universes, so clean up and speed up
            delete h_universe.at(uni).at(k);
        }
    }

    // Loop over the histograms
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        
        // loop over the bins
        for (unsigned int bin = 0; bin < sys_err.at(k).size(); bin ++){
            sys_err.at(k).at(bin) = std::sqrt(sys_err.at(k).at(bin) / universes);
        }
    }

    // Now we have the systematic error computed we should now set the bin error of the CV clone to the std
    // Loop over universes
    for (int uni = 0; uni < universes; uni++){
        
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

        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_universe_2D.at(k)->SetTitle(_util.var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_mcxsec_smear)
            h_universe_2D.at(k)->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
        else if (k == k_xsec_eff)
            h_universe_2D.at(k)->SetTitle(_util.var_labels_eff.at(var).c_str());
        else
            h_universe_2D.at(k)->SetTitle(_util.var_labels_events.at(var).c_str());


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

        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_err->SetTitle(_util.var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_mcxsec_smear || k == k_xsec_mcxsec_shape)
             h_err->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
        else if (k == k_xsec_eff)
            h_err->SetTitle(_util.var_labels_eff.at(var).c_str());
        else
            h_err->SetTitle(_util.var_labels_events.at(var).c_str());

        h_err->SetTitle(" ");
        h_err->SetMarkerColor(kBlack);
        h_err->SetLineStyle(1);
        h_err->SetLineColor(kBlack);

        
        h_err->GetYaxis()->SetTitle("\% Uncertainty");
        h_err->SetMarkerSize(4);
        if (k != k_xsec_dirt && k != k_xsec_ext) h_err->Draw("hist, text00");
        gStyle->SetPaintTextFormat("4.1f");

        FillSysVector(label, var, k, h_err, h_err, h_err);

        gStyle->SetPalette(56);
        gStyle->SetPalette(kBlueGreenYellow);

        if (xsec_types.at(k) != "data_xsec") _util.Draw_ubooneSim(c, 0.40, 0.915, 0.40, 0.915);
        else _util.Draw_Data_POT(c, Data_POT, 0.50, 0.915, 0.50, 0.915);

        // Choose what histograms to save
        if ( (var == k_var_recoX && (k == k_xsec_mcxsec || k == k_xsec_dataxsec || k == k_xsec_bkg || k == k_xsec_sig)) || 
             (var == k_var_trueX && (k == k_xsec_eff || k == k_xsec_mcxsec_smear || k == k_xsec_mcxsec_shape)) ||
             (var == k_var_integrated && (k == k_xsec_mcxsec || k == k_xsec_dataxsec || k == k_xsec_bkg || k == k_xsec_sig || k == k_xsec_eff)) 
            ){
                c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), _util.vars.at(var).c_str(),  _util.run_period, label.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));
            }

        delete c;
        delete leg;
        delete h_universe_2D.at(k);
        delete h_err;

    }
    
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareCVXSec(){

    std::vector<std::string> error_type = {"stat", "sys", "tot"};

    // Loop over the variables
    for (unsigned int var = 0; var < _util.vars.size(); var++){

        // Loop over the error labels
        for (unsigned int err_lab = 0; err_lab < error_type.size(); err_lab++){

            TH1D* h_dataxsec     = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_temp");
            TH1D* h_dataxsec_tot = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_tot_temp"); // For total uncertainty
            TH1D* h_mcxsec       = (TH1D*) cv_hist_vec.at(var).at(k_xsec_mcxsec)  ->Clone("h_mc_xsec_temp");

            TPad *topPad;
            TPad *bottomPad;

            TCanvas * c = new TCanvas("c", "c", 500, 500);
            topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
            bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
            _util.SetTPadOptions(topPad, bottomPad);

            h_dataxsec->SetLineColor(kBlack);
            h_mcxsec  ->SetLineColor(kRed+2);

            // h_dataxsec->GetYaxis()->SetRangeUser(0, 0.5e-39);
            if (_util.vars.at(var) == "integrated") h_dataxsec->GetYaxis()->SetRangeUser(3.5, 10.5);
            else h_dataxsec->GetYaxis()->SetRangeUser(0.0, _util.xsec_scale);

            h_dataxsec->GetYaxis()->SetLabelSize(0.04);
            h_dataxsec->GetYaxis()->SetTitleSize(14);
            h_dataxsec->GetYaxis()->SetTitleFont(44);
            h_dataxsec->GetYaxis()->SetTitleOffset(1.5);
            h_dataxsec->GetXaxis()->SetTitle("");
            h_dataxsec->GetXaxis()->SetLabelSize(0);
            h_dataxsec->SetMarkerStyle(20);
            h_dataxsec->SetMarkerSize(0.5);
            h_dataxsec_tot->SetMarkerStyle(20);
            h_dataxsec_tot->SetMarkerSize(0.5);
            h_dataxsec->SetMinimum(0);
            
            // Rewrite the errors for data to sys
            if (error_type.at(err_lab) == "sys"){
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec->SetBinError(bin+1, 0.01*std::sqrt(v_err.at(k_err_sys).at(var).at(k_xsec_mcxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            // Overwrite error to stat + sys
            else {
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec_tot->SetBinError(bin+1, 0.01*std::sqrt(v_err.at(k_err_sys).at(var).at(k_xsec_mcxsec).at(bin) + v_err.at(k_err_stat).at(var).at(k_xsec_dataxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
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
            if (error_type.at(err_lab) == "stat")     leg->AddEntry(h_dataxsec_tot, "Data (Stat.)", "ep");
            else if (error_type.at(err_lab) == "sys") leg->AddEntry(h_dataxsec_tot, "Data (Sys.)", "ep");
            else                                      leg->AddEntry(h_dataxsec_tot, "Data (Stat.+Sys.)", "ep");
            leg->AddEntry(h_mcxsec_clone,   "MC (Stat.)", "lf");
            leg->Draw();


            bottomPad->cd();
                
            // The percent difference of mc wrt data
            TH1D* h_err = (TH1D *)h_dataxsec->Clone("h_ratio");
            h_err->Add(h_mcxsec, -1);
            h_err->Divide(h_dataxsec);
            h_err->Scale(100);
            h_err->SetLineWidth(2);
            h_err->SetLineColor(kGreen+2);


            bottomPad->SetGridy(kFALSE);
            
            SetRatioOptions(h_err);
            h_err->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
            h_err->SetLineColor(kBlack);
            h_err->GetYaxis()->SetTitleSize(11);
            h_err->GetYaxis()->SetRangeUser(-100, 100);

            // Set the Titles
            h_err->SetTitle(_util.var_labels_xsec.at(var).c_str());
            h_err->GetYaxis()->SetTitle("Data - MC / Data [\%]");
            h_err->GetYaxis()->SetTitleOffset(2.5);
            // h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            if (_util.vars.at(var) == "integrated")  h_err->GetXaxis()->SetLabelSize(0);
            h_err->SetMarkerSize(3.0);
            h_err->Draw("hist,text00");

            // Draw the run period on the plot
            _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

            _util.Draw_Data_POT(c, Data_POT, 0.47, 0.915, 0.47, 0.915);

            
            c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison_%s.pdf", _util.run_period, _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), error_type.at(err_lab).c_str() ));
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

    int mc_xsec_var = k_xsec_mcxsec;

    // Loop over the variables
    for (unsigned int var = 0; var < _util.vars.size(); var++){

        int temp_var = var;

        // For true var, switch comparisons to shape
        if (var == k_var_trueX){
            temp_var = k_var_recoX;
            mc_xsec_var = k_xsec_mcxsec_shape;
        }

        // Loop over the error labels
        for (unsigned int err_lab = 0; err_lab < error_type.size(); err_lab++){

            TH1D* h_dataxsec      = (TH1D*) cv_hist_vec.at(temp_var).at(k_xsec_dataxsec)->Clone("h_data_xsec_temp");
            TH1D* h_dataxsec_tot  = (TH1D*) cv_hist_vec.at(temp_var).at(k_xsec_dataxsec)->Clone("h_data_xsec_tot_temp"); // For total uncertainty
            TH1D* h_mcxsec        = (TH1D*) cv_hist_vec.at(temp_var).at(k_xsec_mcxsec)  ->Clone("h_mc_xsec_temp");
            
            TPad *topPad;
            TPad *bottomPad;

            TCanvas * c = new TCanvas("c", "c", 500, 500);

            h_dataxsec->SetLineColor(kBlack);
            h_mcxsec  ->SetLineColor(kRed+2);

            // h_dataxsec->GetYaxis()->SetRangeUser(0, 0.5e-39);
            if (_util.vars.at(var) == "integrated")
                h_dataxsec->GetYaxis()->SetRangeUser(3.5, 10.5);
            else
                h_dataxsec->GetYaxis()->SetRangeUser(0.0, _util.xsec_scale);

            _util.IncreaseLabelSize(h_dataxsec, c);
            if (_util.vars.at(var) == "integrated")h_dataxsec->GetXaxis()->SetLabelSize(0);
            h_dataxsec->GetYaxis()->SetTitleSize(0.04);
            h_dataxsec->SetTitle(_util.var_labels_xsec.at(var).c_str());
            h_dataxsec->SetMarkerStyle(20);
            h_dataxsec->SetMarkerSize(0.5);
            h_dataxsec_tot->SetMarkerStyle(20);
            h_dataxsec_tot->SetMarkerSize(0.5);
            h_dataxsec->SetMinimum(0);

            // Get the chi-squared
            double chi, pval;
            int ndof;
            int index = 0;
            if (error_type.at(err_lab) == "stat"){
                index = k_err_stat;
            }
            else if (error_type.at(err_lab) == "sys"){
                index = k_err_sys;
            }
            else{
                index = k_err_tot;
            }
            
            TH2D *h_cov_temp = (TH2D*)h_cov_v.at(var).at(mc_xsec_var).at(index)->Clone();
        
            TH2D *h_cov_temp_data_stat = (TH2D*)h_cov_v.at(var).at(k_xsec_dataxsec).at(k_err_stat)->Clone();;

            // Need to convert the cov units for stat
            if (var == k_var_trueX){
                _util.ConvertCovarianceUnits(h_cov_temp_data_stat, 
                        cv_hist_vec.at(var).at(k_xsec_dataxsec), 
                        h_mcxsec);
            }

            if (index == k_err_stat)
                h_cov_temp->Add(h_cov_temp_data_stat);


            // Calculate the chi-squared before plotting
            _util.CalcChiSquared(h_mcxsec, h_dataxsec_tot, h_cov_temp , chi, ndof, pval);  

            // Convert to data units
            _util.ConvertCovarianceUnits(h_cov_temp, 
                        cv_hist_vec.at(var).at(mc_xsec_var), 
                        cv_hist_vec.at(temp_var).at(k_xsec_dataxsec));

            // Rewrite the errors for data to sys
            if (error_type.at(err_lab) == "sys"){
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec->SetBinError(bin+1, std::sqrt(h_cov_temp->GetBinContent(bin+1, bin+1)) );

                    std::cout <<"Sys: "<< 100*h_dataxsec->GetBinError(bin+1) / h_dataxsec->GetBinContent(bin+1) << std::endl;
                }
            }
            else if (error_type.at(err_lab) == "stat"){
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    // h_dataxsec->SetBinError(bin+1, std::sqrt(h_cov_temp->GetBinContent(bin+1, bin+1)));

                    std::cout <<"Stat: " <<100*h_dataxsec->GetBinError(bin+1) / h_dataxsec->GetBinContent(bin+1) << std::endl;
                }
            }
            // Overwrite error to stat + sys
            else {
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec_tot->SetBinError(bin+1, std::sqrt(h_cov_temp->GetBinContent(bin+1, bin+1)));
                    // h_dataxsec_tot->SetBinError(bin+1, 0.257 * h_dataxsec->GetBinContent(bin+1)  );

                    std::cout << "Tot: "<< 100*h_dataxsec_tot->GetBinError(bin+1) / h_dataxsec_tot->GetBinContent(bin+1) << std::endl;
                }

            }

            if (_util.zoom && std::string(_util.xsec_var) == "elec_cang"){
                h_dataxsec->GetXaxis()->SetRangeUser(0.6, 1.0);
            }

            h_dataxsec->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
            if (var == k_var_integrated){
                h_dataxsec->GetXaxis()->SetLabelSize(0);
                h_dataxsec->SetTitle("");
            }

            c->SetLeftMargin(0.2);
            c->SetBottomMargin(0.15);
            h_dataxsec->GetYaxis()->SetTitleOffset(1.8);
            h_dataxsec->Draw("E,X0");
            if (error_type.at(err_lab) == "tot"){
                h_dataxsec_tot->Draw("E1,same,X0");
                h_dataxsec->Draw("E1,same,X0");

            } 

            h_mcxsec->Draw("hist,same");
            
            TH1D* h_mcxsec_clone = (TH1D *)h_mcxsec->Clone("h_mc_clone");
            h_mcxsec_clone->SetFillColorAlpha(0, 0);
            h_mcxsec_clone->Draw("hist,same");

            // redraw the data so the data is on top of everything
            h_dataxsec->Draw("E,X0,same");
            if (error_type.at(err_lab) == "tot"){
                h_dataxsec_tot->Draw("E1,same,X0,same");
                h_dataxsec->Draw("E1,same,X0,same");

            }

                      

            TLegend *leg = new TLegend(0.31, 0.6, 0.66, 0.85);
            gStyle->SetLegendTextSize(0.03);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if (error_type.at(err_lab) == "stat")     leg->AddEntry(h_dataxsec_tot, "Data (Stat.)", "ep");
            else if (error_type.at(err_lab) == "sys") leg->AddEntry(h_dataxsec_tot, "Data (Sys.)", "ep");
            else                                      leg->AddEntry(h_dataxsec_tot, "Data (Stat. + Sys.)", "ep");
            leg->AddEntry(h_mcxsec_clone,   Form("GENIE v3.0.6 (#muB tune) #chi^{2}/N_{dof} = %2.1f/%i", chi, ndof), "l");
            leg->Draw();

            // Draw the run period on the plot
            // _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

            _util.Draw_Data_POT(c, Data_POT, 0.52, 0.92, 0.52, 0.92);

            if (std::string(_util.xsec_smear_mode) == "er"){
                if (_util.zoom)
                    c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison_%s_no_ratio_zoom.pdf", _util.run_period, _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), error_type.at(err_lab).c_str() ));
                else
                    c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison_%s_no_ratio.pdf", _util.run_period, _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), error_type.at(err_lab).c_str() ));
            }
            delete c;
            delete h_dataxsec;
            delete h_mcxsec;
            delete h_cov_temp;

        
        }

    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialsePlotCV(){

    // Get the CV histograms. These should stay constant througout the code

    cv_hist_vec.resize(_util.vars.size());
    
    for (unsigned int var = 0; var < _util.vars.size(); var++){
        cv_hist_vec.at(var).resize(xsec_types.size());
    }


    // Loop over the _util.vars
    for (unsigned int var = 0; var < _util.vars.size(); var++){
        
        // Loop over the typrs
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            _util.GetHist(f_nuexsec, cv_hist_vec.at(var).at(k), Form( "CV/%s/h_run%s_CV_0_%s_%s", _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

            if (cv_hist_vec.at(var).at(k) == NULL) std::cout << "Failed to get the histogram!" << std::endl;

            if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
                cv_hist_vec.at(var).at(k)->Scale(1.0, "width");

            // Customise
            cv_hist_vec.at(var).at(k)->SetLineWidth(2);
            cv_hist_vec.at(var).at(k)->SetLineColor(kBlack);

            // Set the Titles
            if (k == k_xsec_mcxsec)
                cv_hist_vec.at(var).at(k)->SetTitle(_util.var_labels_xsec.at(var).c_str());
            else if (k == k_xsec_mcxsec_smear)
                cv_hist_vec.at(var).at(k)->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
            else if (k == k_xsec_dataxsec)
                cv_hist_vec.at(var).at(k)->SetTitle(_util.var_labels_events.at(var).c_str());
            else if (k == k_xsec_eff)
                cv_hist_vec.at(var).at(k)->SetTitle(_util.var_labels_eff.at(var).c_str());
            else
                cv_hist_vec.at(var).at(k)->SetTitle(_util.var_labels_events.at(var).c_str());

        }

        // Create the CV directory and draw the CV
        _util.CreateDirectory("/Systematics/CV/" + _util.vars.at(var) + "/");
    }
    
    TCanvas *c_cv;
    
    // Loop over the _util.vars
    for (unsigned int var = 0; var < _util.vars.size(); var++){
        
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

            c_cv->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_%s.pdf", _util.run_period, _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

            delete c_cv;
            delete leg;
        }
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareVariationXSec(std::string label, int var, std::string label_pretty){

    
    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + _util.vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe;
    
    // Resize to the number of universes
    h_universe.resize(2);
    h_universe.at(k_up).resize(xsec_types.size());
    h_universe.at(k_dn).resize(xsec_types.size());

    // Now get the histograms
    std::string label_up = label + "up";
    std::string label_dn = label + "dn";
    
    TH1D* htemp;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), _util.vars.at(var).c_str(), _util.run_period, label_up.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

        h_universe.at(k_up).at(k) = (TH1D*)htemp->Clone();

        if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
            h_universe.at(k_up).at(k)->Scale(1.0, "width");

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        if (k == k_xsec_mcxsec || k == k_xsec_mcxsec_smear) h_universe.at(k_up).at(k)->SetLineStyle(7);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), _util.vars.at(var).c_str(), _util.run_period, label_dn.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

        h_universe.at(k_dn).at(k) = (TH1D*)htemp->Clone();

        if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
            h_universe.at(k_dn).at(k)->Scale(1.0, "width");

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
        if (k == k_xsec_mcxsec || k == k_xsec_mcxsec_smear) h_universe.at(k_dn).at(k)->SetLineStyle(7);
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
    if (_util.vars.at(var) != "integrated"){
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

    // Set the Titles
    if (var == k_xsec_mcxsec)
        h_err_up->SetTitle(_util.var_labels_xsec.at(var).c_str());
    else if (var == k_xsec_dataxsec)
        h_err_up->SetTitle(_util.var_labels_events.at(var).c_str());
    else if (var == k_xsec_eff)
        h_err_up->SetTitle(_util.var_labels_eff.at(var).c_str());
    else
        h_err_up->SetTitle(_util.var_labels_events.at(var).c_str());

    // h_err_up->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
    
    h_err_up->Draw("hist,same");
    h_err_dn->Draw("hist,same");
    h_err->Draw("hist,same");

    c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_data_mc_comparison.pdf", _util.run_period, label.c_str(), _util.vars.at(var).c_str(), _util.run_period, label.c_str(), _util.vars.at(var).c_str() ));

    // Clear memory
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        delete h_universe.at(k_up).at(k);
        delete h_universe.at(k_dn).at(k);
    }

    delete c;
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CalcMatrices(std::string label, int var, std::vector<std::vector<TH1D*>> h_universe, int _type, TH1D* h_CV ){


    int n_bins = cv_hist_vec.at(var).at(0)->GetNbinsX();

    TH2D* cov; 

    if (var == k_var_integrated){
        cov  = new TH2D("h_cov",   ";Bin i; Bin j",1, 1, 2, 1, 1, 2 ); // Create the covariance matrix
    }
    else { 
        cov  = new TH2D("h_cov",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the covariance matrix
    }


    // Loop over universes
    // std::cout << "Universes: " << h_universe.size() << std::endl;
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over the rows
        for (int row = 1; row < h_CV->GetNbinsX()+1; row++){
            
            double uni_row = h_universe.at(uni).at(_type)->GetBinContent(row);
            double cv_row  = h_CV->GetBinContent(row);

            // Loop over the columns
            for (int col = 1; col < h_CV->GetNbinsX()+1; col++){

                double uni_col = h_universe.at(uni).at(_type)->GetBinContent(col);
                double cv_col  = h_CV->GetBinContent(col);
                
                double c = (uni_row - cv_row) * (uni_col - cv_col);

                if (uni != h_universe.size()-1)    cov->SetBinContent(row, col, cov->GetBinContent(row, col) + c ); // Fill with variance 
                else cov->SetBinContent(row, col, (cov->GetBinContent(row, col) + c) / h_universe.size());       // Fill with variance and divide by nuni
            
            }

        }
            
    }

    gStyle->SetPalette(kBlueGreenYellow);
    
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
        label == "NormNCCOH"        ||
        label == "xsr_scc_Fv3"      ||
        label == "xsr_scc_Fa3"){
            h_cov_v.at(var).at(_type).at(k_err_genie_uni)->Add(cov);
            h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
            h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    // Beamline
    else if (label == "Horn_curr" ||
            label == "Horn1_x" ||
            label == "Horn1_y" ||
            label == "Beam_spot" ||
            label == "Horn2_x" ||
            label == "Horn2_y" ||
            label == "Horn_Water" ||
            label == "Beam_shift_x" ||
            label == "Beam_shift_y" ||
            label == "Target_z" ){
            h_cov_v.at(var).at(_type).at(k_err_beamline)->Add(cov);
            h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
            h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);
    }
    // Detector Variation
    else if (
            label == "LYRayleigh" ||
            label == "LYDown"     ||
            label == "LYAttenuation" ||
            label == "SCE" ||
            label == "Recomb2" ||
            label == "WireModX" ||
            label == "WireModYZ" ||
            label == "WireModThetaXZ" ||
            label == "WireModThetaYZ_withSigmaSplines" ||
            label == "WireModThetaYZ_withoutSigmaSplines" ||
            label == "WireModdEdX" ){
            
            // Convert the Covariance Matrix-- switching from detvar cv deviations to CV deviation
            _util.ConvertCovarianceUnits(cov, h_CV, cv_hist_vec.at(var).at(_type));
           
            h_cov_v.at(var).at(_type).at(k_err_detvar)->Add(cov);
            h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
            h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);
    }
    else if (label == "weightsGenie"){
        h_cov_v.at(var).at(_type).at(k_err_genie_multi)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    else if (label == "weightsReint"){
        h_cov_v.at(var).at(_type).at(k_err_reint)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    else if (label == "weightsFlux"){
        h_cov_v.at(var).at(_type).at(k_err_hp)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    else if (label == "Dirt"){
        h_cov_v.at(var).at(_type).at(k_err_dirt)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    else if (label == "POT"){
        h_cov_v.at(var).at(_type).at(k_err_pot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    else if (label == "MCStats"){
        h_cov_v.at(var).at(_type).at(k_err_mcstats)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    else if (label == "EXT"){
        h_cov_v.at(var).at(_type).at(k_err_ext)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_tot)->Add(cov);
        h_cov_v.at(var).at(_type).at(k_err_sys)->Add(cov);

    }
    else {
        std::cout << "Unknown variation specified: " << label << std::endl;
        delete cov;
        return;
    }

    delete cov;

}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillSysVector(std::string variation, int var, int type, TH1D *h_up, TH1D *h_dn, TH1D* h_CV){

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
        variation == "NormNCCOH"        ||
        variation == "xsr_scc_Fv3"      ||
        variation == "xsr_scc_Fa3"){

        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Error is given as the square uncertainty in percent

            // Note only RPA is a two univerese unisim, the others are single variations
            if (variation == "RPA"){
                av_err += h_up->GetBinContent(bin+1) * h_up->GetBinContent(bin+1);
                av_err += h_dn->GetBinContent(bin+1) * h_dn->GetBinContent(bin+1);
                av_err = std::sqrt(av_err / 2.0 );
            }
            else {
                av_err += std::sqrt(h_up->GetBinContent(bin+1) * h_up->GetBinContent(bin+1));
            }

            av_err = 100 * av_err / h_CV->GetBinContent(bin+1);

            v_err.at(k_err_genie_uni).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin)       += av_err*av_err;
            
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
        variation == "Target_z" ){

        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Error is given as the square uncertainty in percent
            
            // Singular variation
            if (variation == "Decay_pipe_Bfield"){
                av_err += std::sqrt(h_up->GetBinContent(bin+1) * h_up->GetBinContent(bin+1));
            }
            // Two sided variation
            else {
                av_err += h_up->GetBinContent(bin+1) * h_up->GetBinContent(bin+1);
                av_err += h_dn->GetBinContent(bin+1) * h_dn->GetBinContent(bin+1);
                av_err = std::sqrt(av_err / 2.0 );
            }

            av_err = 100 * av_err / h_CV->GetBinContent(bin+1);

            v_err.at(k_err_beamline).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin)      += av_err*av_err;
            
        }
    }
    // This is a detector variation
    else if (
        variation == "LYRayleigh" ||
        variation == "LYDown"     ||
        variation == "LYAttenuation" ||
        variation == "SCE" ||
        variation == "Recomb2" ||
        variation == "WireModX" ||
        variation == "WireModYZ" ||
        variation == "WireModThetaXZ" ||
        variation == "WireModThetaYZ_withSigmaSplines" ||
        variation == "WireModThetaYZ_withoutSigmaSplines" ||
        variation == "WireModdEdX"
        ){

        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            
            // Error is given as the square uncertainty in percent
            av_err += std::sqrt(h_up->GetBinContent(bin+1) * h_up->GetBinContent(bin+1));
            av_err = 100 * av_err / h_CV->GetBinContent(bin+1);
            v_err.at(k_err_detvar).at(var).at(type).at(bin)   += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin)      += av_err*av_err;
            
        }
    }
    else if (variation == "weightsGenie"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Get the error in each bin, then add the square
            av_err += std::abs(h_up->GetBinContent(bin+1));
            v_err.at(k_err_genie_multi).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "weightsReint"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Error is given as the square uncertainty in percent
            av_err += std::abs(h_up->GetBinContent(bin+1));
            v_err.at(k_err_reint).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "weightsFlux"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Error is given as the square uncertainty in percent
            av_err += std::abs(h_up->GetBinContent(bin+1));
            v_err.at(k_err_hp).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "MCStats"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Error is given as the square uncertainty in percent
            av_err += std::abs(h_up->GetBinContent(bin+1));
            v_err.at(k_err_mcstats).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "Dirt"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Error is given as the square uncertainty in percent
            av_err += std::abs(h_up->GetBinContent(bin+1));
            av_err = 100 * av_err / h_CV->GetBinContent(bin+1);
            v_err.at(k_err_dirt).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin)  += av_err*av_err;
            
        }

    }
    else if (variation == "POT"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
           // Error is given as the square uncertainty in percent
            av_err += std::abs(h_up->GetBinContent(bin+1));
            av_err = 100 * av_err / h_CV->GetBinContent(bin+1);
            v_err.at(k_err_pot).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else if (variation == "EXT"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double av_err = 0;
            // Error is given as the square uncertainty in percent
            av_err += std::abs(h_up->GetBinContent(bin+1));
            av_err = 100 * av_err / h_CV->GetBinContent(bin+1);
            v_err.at(k_err_ext).at(var).at(type).at(bin) += av_err*av_err;
            v_err.at(k_err_sys).at(var).at(type).at(bin) += av_err*av_err;
            
        }

    }
    else {
        std::cout << "Unknown variation specified: " << variation << std::endl;
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillStatVector(){

    // Loop over the differential variables
    for (unsigned int var = 0; var < _util.vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            
            // Get the uncertainty in each bin
            for (int bin = 0; bin < cv_hist_vec.at(var).at(type)->GetNbinsX(); bin++){
            
                double stat_err = 100 * cv_hist_vec.at(var).at(type)->GetBinError(bin+1) / cv_hist_vec.at(var).at(type)->GetBinContent(bin+1);
                v_err.at(k_err_stat).at(var).at(type).at(bin) += stat_err*stat_err;

                // To avoid double counting the N - B uncertainty, set the MC xsec stat error to zero.
                if ((var == k_var_recoX || var == k_var_trueX) && (type == k_xsec_mcxsec || type == k_xsec_mcxsec_shape)){
                    cv_hist_vec.at(var).at(type)->SetBinError(bin+1,0);
                    v_err.at(k_err_stat).at(var).at(type).at(bin) = 0.0;
                }

            }
        
        }
    }

    // Lets also fill the diagonals of the statistical covariance matrix

     // Loop over the differential variables
    for (unsigned int var = 0; var < _util.vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){

            // loop over rows
            for (int row = 1; row < h_cov_v.at(var).at(type).at(k_err_stat)->GetNbinsX()+1; row++) {
                
                // Loop over columns
                for (int col = 1; col < h_cov_v.at(var).at(type).at(k_err_stat)->GetNbinsY()+1; col++) {
                    
                    // Only set the bin content of the diagonals
                    double bin_diag = cv_hist_vec.at(var).at(type)->GetBinContent(row); // We use the MC value to fill the cov matrix, but use the data stat err for now. 

                    // 0.01 converts each percentage back to a number. We multiply this by the cv to get the deviate
                    if (row == col) h_cov_v.at(var).at(type).at(k_err_stat)->SetBinContent(row, col, 0.01*0.01*v_err.at(k_err_stat).at(var).at(type).at(row-1)*bin_diag*bin_diag);  
                }
            }

            // Add the created stat diagonal cov matrix to the total
            h_cov_v.at(var).at(type).at(k_err_tot)->Add(h_cov_v.at(var).at(type).at(k_err_stat));
        }
    }

    // ---

    // Add the data stat to the total mc cov matrix
    // Convert the Covariance Matrix-- switching from detvar Data CV deviations to MC CV deviation
    _util.ConvertCovarianceUnits(h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(k_err_stat), 
                           cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec), 
                           cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec));
    
    h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_tot)->Add(h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(k_err_stat));

    // Put it back
    // Convert the Covariance Matrix-- switching from detvar MC CV deviations to data CV deviation
    _util.ConvertCovarianceUnits(h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(k_err_stat), 
                           cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec), 
                           cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec));

    // ---

    // Add the MC stat err tot the total MC sys error
    h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_sys)->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_stat));

    // ---

    // Add the data stat errors to the shape channel total covariance matrix
    // Convert the Covariance Matrix-- switching from Data CV deviations to MC CV deviation
    _util.ConvertCovarianceUnits(h_cov_v.at(k_var_trueX).at(k_xsec_dataxsec).at(k_err_stat), 
                           cv_hist_vec.at(k_var_trueX).at(k_xsec_dataxsec), 
                           cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape));

    h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot)->Add(h_cov_v.at(k_var_trueX).at(k_xsec_dataxsec).at(k_err_stat));
    
    // put it back
    // Convert the Covariance Matrix-- switching from mc shape CV deviations to data deviation
    _util.ConvertCovarianceUnits(h_cov_v.at(k_var_trueX).at(k_xsec_dataxsec).at(k_err_stat), 
                           cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape),
                           cv_hist_vec.at(k_var_trueX).at(k_xsec_dataxsec));

}
// -----------------------------------------------------------------------------
void SystematicsHelper::AddSmearCovMatrix(){

    // Only do this in er and wiener modes
    if (std::string(_util.xsec_smear_mode) == "er" || std::string(_util.xsec_smear_mode) == "wiener"){

        // Loop over the error matrices and add the smear uncertainty in
        for (unsigned int cov = 0; cov < h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).size(); cov++){

            // Skip the stat err
            if (cov == k_err_stat)
                continue;

            // MC event rate
            h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(cov)->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_smear).at(cov));

            // --

            // Convert the Covariance Matrix-- switching from MC CV deviations to Data CV deviation
            _util.ConvertCovarianceUnits(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_smear).at(cov),
                               cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec),
                               cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec));

            // Data event rate
            h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(cov)->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_smear).at(cov));

            // Put it back
            // Convert the Covariance Matrix-- switching from MC CV deviations to Data CV deviation
            _util.ConvertCovarianceUnits(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_smear).at(cov),
                               cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec),
                               cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec));
        }
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillPOTCountingVector(){

    // Loop over the differential variables
    for (unsigned int var = 0; var < _util.vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            
            // Get the uncertainty in each bin
            for (int bin = 0; bin < cv_hist_vec.at(var).at(type)->GetNbinsX(); bin++){
            
                v_err.at(k_err_pot).at(var).at(type).at(bin) += 2.0*2.0;
                v_err.at(k_err_sys).at(var).at(type).at(bin) += 2.0*2.0;

            }
        
        }
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PrintUncertaintySummary(){

    // Loop over the variables
    for (unsigned int var = 0; var < _util.vars.size(); var++ ){
        // if (var == k_var_recoX) continue; // Skip the true var which doesnt make much sense

        std::cout <<"----------------------------------------------" << std::endl;
        std::cout <<"Differential Variable: " << _util.vars.at(var) <<"\n"<< std::endl;
        
        // Loop over the bins
        for (unsigned int bin = 0; bin < v_err.front().at(var).at(k_xsec_mcxsec).size(); bin++ ){


            double bin_diag_mc = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(bin+1);
            double bin_diag_mc_shape = cv_hist_vec.at(var).at(k_xsec_mcxsec_shape)->GetBinContent(bin+1);
            double bin_diag_data = cv_hist_vec.at(var).at(k_xsec_dataxsec)->GetBinContent(bin+1);
            double bin_diag = bin_diag_mc_shape;

            int _type = k_xsec_mcxsec_shape;

            if (var == k_var_integrated){
                _type = k_xsec_dataxsec;
                bin_diag = bin_diag_data;
            }

            if (var == k_var_recoX){
                _type = k_xsec_mcxsec;
                bin_diag = bin_diag_mc;
            }
            
            std::cout << "Bin: " << bin+1 << " GENIE Unisim:       " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_genie_uni)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " GENIE Multisim:     " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_genie_multi)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Beamline:           " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_beamline)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Hadron Prod.:       " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_hp)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Geant Rein.:        " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_reint)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Detector:           " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_detvar)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Dirt:               " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_dirt)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " POT Counting:       " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_pot)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " MC Stats:           " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_mcstats)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " EXT Norm:           " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_ext)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            
            // if ( (std::string(_util.xsec_smear_mode) == "er" || std::string(_util.xsec_smear_mode) == "wiener") && var == k_var_recoX)
            //     std::cout << "Bin: " << bin+1 << " Smearing:           " << std::sqrt(v_err.at(k_err_sys).at(k_var_trueX).at(k_xsec_mcxsec_smear).at(bin)) << " \%"<< std::endl;

            std::cout << std::endl;

            // Statistical Uncertainties
            std::cout << "Bin: " << bin+1 << " Tot Data X-Sec Stat:                 " << 100 * std::sqrt(h_cov_v.at(var).at(k_xsec_dataxsec).at(k_err_stat)->GetBinContent(bin+1, bin+1)) / bin_diag_data << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Stat:                   " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_mcstats)->GetBinContent(bin+1, bin+1)) / bin_diag_mc    << " \%"<< std::endl;
            
            // Systematic Uncertainties
            if (var == k_var_integrated)
                std::cout << "Bin: " << bin+1 << " Tot Data X-Sec Sys:                  " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_sys)->GetBinContent(bin+1, bin+1)) / bin_diag  << " \%"<< std::endl;
            else
                std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Sys:                    " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_sys)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            
            // Total Uncertanities
            if ( (std::string(_util.xsec_smear_mode) == "er" || std::string(_util.xsec_smear_mode) == "wiener") && var == k_var_recoX){
                std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Uncertainty:            " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_tot)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            }
            else {
                if (var == k_var_integrated)
                    std::cout << "Bin: " << bin+1 << " Tot Data X-Sec Uncertainty:          " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_tot)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
                else
                    std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Uncertainty:            " << 100 * std::sqrt(h_cov_v.at(var).at(_type).at(k_err_tot)->GetBinContent(bin+1, bin+1)) / bin_diag << " \%"<< std::endl;
            }
            std::cout <<"\n--" << std::endl;       
        }
    }
    std::cout <<"----------------------------------------------" << std::endl;

    // Print the sqrt of the diagonals of the covariance matrix -- should be equivalent
    bool print_debug = false;
    if (print_debug){
        // loop over rows
        for (int i = 1; i < h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_tot)->GetNbinsX()+1; i++) {
            
                // Only set the bin content of the diagonals
                double bin_diag = cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec)->GetBinContent(i);
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_genie_uni)->GetBinContent(i, i)) / bin_diag   << " genie_uni" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_genie_multi)->GetBinContent(i, i)) / bin_diag << " genie_multi" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_hp)->GetBinContent(i, i)) / bin_diag          << " hp" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_beamline)->GetBinContent(i, i)) / bin_diag    << " beamline" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_dirt)->GetBinContent(i, i)) / bin_diag        << " dirt" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_pot)->GetBinContent(i, i)) / bin_diag         << " pot" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_reint)->GetBinContent(i, i)) / bin_diag       << " reint" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_detvar)->GetBinContent(i, i)) / bin_diag      << " detvar" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_ext)->GetBinContent(i, i)) / bin_diag         << " EXT" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_mcstats)->GetBinContent(i, i)) / bin_diag     << " mcstats" << std::endl;
                
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(k_err_stat)->GetBinContent(i, i)) / bin_diag << " data stat" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_sys)->GetBinContent(i, i)) / bin_diag << " sys" << std::endl;
                std::cout << "Bin" << i << ": "<< 100 * std::sqrt(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_tot)->GetBinContent(i, i)) / bin_diag << " tot" << std::endl;
                std::cout << std::endl;
        }
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseUncertaintyVectors(){

    // differential variable, type, bin error
    v_err.resize(k_ERR_MAX);

    // Loop over the systematic error types
    for (unsigned int err = 0; err < v_err.size(); err++){
        v_err.at(err).resize(_util.vars.size());
    }

    // Loop over the systematic error types
    for (unsigned int err = 0; err < v_err.size(); err++){

        // Loop over the _util.vars
        for (unsigned int var = 0; var < _util.vars.size(); var++ ){
            v_err.at(err).at(var).resize(xsec_types.size());
        }
    }

    // Loop over the systematic error types
    for (unsigned int err = 0; err < v_err.size(); err++){
        
        // Loop over the _util.vars
        for (unsigned int var = 0; var < _util.vars.size(); var++ ){
            
            // Loop over the types
            for (unsigned int type = 0; type < xsec_types.size(); type++ ){
                v_err.at(err).at(var).at(type).resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            }
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
    if (std::string(_util.xsec_var) == "elec_E"){
            cov->GetZaxis()->SetTitle("Covariance [cm^{4}/GeV^{2}]");
    }
    if (std::string(_util.xsec_var) == "elec_ang"){
        cov->GetZaxis()->SetTitle("Covariance [cm^{4}/deg^{2}]");
    }
    if (std::string(_util.xsec_var) == "elec_cang"){
        cov->GetZaxis()->SetTitle("Covariance [cm^{4}]");
        cov->GetXaxis()->SetNdivisions(cov->GetNbinsX(), 0, 0, kFALSE);
        cov->GetYaxis()->SetNdivisions(cov->GetNbinsY(), 0, 0, kFALSE);
    }
    gStyle->SetPaintTextFormat("0.3f");
    gStyle->SetPaintTextFormat("1.1g");
    cov->SetMarkerColor(kRed+1);
    cov->GetZaxis()->SetMaxDigits(2);
    cov->GetZaxis()->SetTitleOffset(1.45);
    cov->SetTitle("Covariance Matrix");
    cov->Draw("colz,text00");
    _util.IncreaseLabelSize(cov, c);
    cov->SetMarkerSize(1.3);
    _util.Draw_ubooneSim(c, 0.30, 0.915, 0.30, 0.915);
    c->Print(print_name.c_str());
    delete c;

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SaveCorMatrix(TH2D* cov, TH1D* h_CV, std::string print_name){

    gStyle->SetPalette(kBlueGreenYellow);

    TH2D* cor = (TH2D*)cov->Clone();
    _util.CalcCorrelation(h_CV, cov, cor);

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    cor->GetXaxis()->CenterLabels(kTRUE);
    cor->GetYaxis()->CenterLabels(kTRUE);
    cor->GetZaxis()->CenterTitle();
    cor->GetZaxis()->SetTitleOffset(1.45);
    cor->SetMarkerColor(kRed+1);
    gStyle->SetPaintTextFormat("0.3f");
    cor->SetTitle("Correlation Matrix");
    cor->GetZaxis()->SetTitle("Correlation");
    if (std::string(_util.xsec_var) == "elec_cang"){
        cor->GetXaxis()->SetNdivisions(cor->GetNbinsX(), 0, 0, kFALSE);
        cor->GetYaxis()->SetNdivisions(cor->GetNbinsY(), 0, 0, kFALSE);
    }
    cor->GetZaxis()->SetRangeUser(-1, 1);
    cor->Draw("colz, text00");
    _util.IncreaseLabelSize(cor, c);
    cor->SetMarkerSize(1.3);
    _util.Draw_ubooneSim(c, 0.30, 0.915, 0.30, 0.915);
    c->Print(print_name.c_str());
    delete c;
    delete cor;

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SaveFracCovMatrix(TH2D* cov, TH1D* h_CV, std::string print_name){

    gStyle->SetPalette(kBlueGreenYellow);

    TH2D* frac_cov= (TH2D*)cov->Clone();
    _util.CalcCFracCovariance(h_CV, frac_cov);

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    frac_cov->GetXaxis()->CenterLabels(kTRUE);
    frac_cov->GetYaxis()->CenterLabels(kTRUE);
    frac_cov->GetZaxis()->CenterTitle();
    frac_cov->GetZaxis()->SetTitleOffset(1.45);
    frac_cov->SetTitle("Fractional Covariance Matrix");
    frac_cov->GetZaxis()->SetTitle("Frac. Covariance");
    frac_cov->SetMarkerColor(kRed+1);
    gStyle->SetPaintTextFormat("0.3f");
    if (std::string(_util.xsec_var) == "elec_cang"){
        frac_cov->GetXaxis()->SetNdivisions(frac_cov->GetNbinsX(), 0, 0, kFALSE);
        frac_cov->GetYaxis()->SetNdivisions(frac_cov->GetNbinsY(), 0, 0, kFALSE);
    }
    gStyle->SetPaintTextFormat("1.1g");
    frac_cov->Draw("colz,text00");
    frac_cov->GetZaxis()->SetMaxDigits(2);
    _util.IncreaseLabelSize(frac_cov, c);
    frac_cov->SetMarkerSize(1.3);
    _util.Draw_ubooneSim(c, 0.30, 0.915, 0.30, 0.915);
    c->Print(print_name.c_str());
    delete c;
    delete frac_cov;

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotTotUnisim(std::string unisim_type){

    std::vector<std::vector<std::vector<double>>> v_unisim;
    bool is_detvar = false;

    std::vector<std::string> unisim_names;

    // Use the beamline errors
    if (unisim_type == "Beamline")   {
        v_unisim = v_err.at(k_err_beamline);
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
                    "Target_z"
                };
    }
    else if (unisim_type == "Genie_Unisim") {
        v_unisim = v_err.at(k_err_genie_uni);

        unisim_names = {
                    "RPA",
                    "CCMEC",
                    "AxFFCCQE",
                    "VecFFCCQE",
                    "DecayAngMEC",
                    "ThetaDelta2Npi",
                    "ThetaDelta2NRad",
                    // "RPA_CCQE_Reduced",
                    "NormCCCOH",
                    "NormNCCOH",
                    "xsr_scc_Fv3",
                    "xsr_scc_Fa3"
                };
    }
    else if (unisim_type == "DetVar") {
        v_unisim = v_err.at(k_err_detvar);
        unisim_names = var_string;
        is_detvar = true;
    }
    else {
        std::cout << "Unknown unisim type specified" << std::endl;
        return;
    }

    // Get all the universes so we can draw them on
    std::vector<std::vector<std::vector<std::vector<TH1D*>>>> h_universe; // var -- label -- Up/Dn -- Type
    
    h_universe.resize(_util.vars.size());
    
    for (unsigned int var = 0; var < h_universe.size(); var++){
        h_universe.at(var).resize(unisim_names.size());
    }
    
    // Loop over the _util.vars
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
    std::vector<bool> single_var;

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

            // Skip the beamline variations we dont want
            if (is_detvar && (unisim_names.at(label) == "CV" || unisim_names.at(label) == "BNB_Diffusion") )
                continue;


            // Check if its just a single on/off type variation
            // This case we dont want to plot the up/dn, but just once
            if (label_up == label_dn){
                single_var.push_back(true);
            }
            else
                single_var.push_back(false);

            // Get the histograms and customise a bit
            for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

                TH1D* htemp;
                
                // Beamline have CV in the name because of the way the cross section helper works
                if (!is_detvar){
                    _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), _util.vars.at(var).c_str(), _util.run_period, label_up.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));
                    
                    if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
                        htemp->Scale(1.0, "width");
                }
                else {
                    _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_CV_0_%s_%s", label_up.c_str(), _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

                    if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
                        htemp->Scale(1.0, "width");
                
                    // Also need to Scale the beamline histograms to the right POT
                    double scale_fact = POT_v.at(k_CV) / POT_v.at(label);
                    
                    // Scale the histograms, but only in the case of certain variables
                    if (xsec_types.at(k) == "sel" || xsec_types.at(k) == "gen" || xsec_types.at(k) == "bkg" || xsec_types.at(k) == "sig")
                        htemp->Scale(scale_fact);
                }

                h_universe.at(var).at(label).at(k_up).at(k) = (TH1D*) htemp->Clone(Form("h_clone_%s_%s_up", _util.vars.at(var).c_str(), labels_up_v.at(label).c_str()));

                // Customise
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineWidth(2);
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineStyle(0);
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineColor(kGreen+2);
                
                // Beamline have CV in the name because of the way the cross section helper works
                if (!is_detvar){
                    _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), _util.vars.at(var).c_str(), _util.run_period, label_dn.c_str(), _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

                    if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
                        htemp->Scale(1.0, "width");
                }
                else {
                    _util.GetHist(f_nuexsec, htemp, Form( "%s/%s/h_run%s_CV_0_%s_%s", label_dn.c_str(), _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(k).c_str()));

                    if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && k != k_xsec_eff ))
                        htemp->Scale(1.0, "width");
                    
                    // Also need to Scale the beamline histograms to the right POT
                    double scale_fact = POT_v.at(k_CV) / POT_v.at(label);
                    
                    // Scale the histograms, but only in the case of MC variables
                    if (xsec_types.at(k) != "ext" && xsec_types.at(k) != "dirt" && xsec_types.at(k) != "data") htemp->Scale(scale_fact);
                }

                h_universe.at(var).at(label).at(k_dn).at(k) = (TH1D*) htemp->Clone(Form("h_clone_%s_%s_up", _util.vars.at(var).c_str(), labels_dn_v.at(label).c_str()));

                // Customise
                h_universe.at(var).at(label).at(k_dn).at(k)->SetLineWidth(2);
                h_universe.at(var).at(label).at(k_dn).at(k)->SetLineStyle(1);
                h_universe.at(var).at(label).at(k_dn).at(k)->SetLineColor(kRed+2);
                h_universe.at(var).at(label).at(k_up).at(k)->SetLineColor(kGreen+2);
            }

           // In the case of a true variation, do the CV background subtraction 
            if (var == k_var_trueX && is_detvar){
                
                TH1D* h_bkg_CV, *h_ext_CV, *h_dirt_CV;
                
                // Load in the backgrounds from file
                TH1D* h_temp = (TH1D*)f_nuexsec->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_bkg", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
                h_bkg_CV = (TH1D*)h_temp->Clone();

                h_temp = (TH1D*)f_nuexsec->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_ext", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
                h_ext_CV = (TH1D*) h_temp->Clone();
                h_ext_CV->Scale(_util.ext_scale_factor /  (_util.config_v.at(_util.k_Run1_Data_POT)  / POT_v.at(k_CV)) );
                h_bkg_CV->Add(h_ext_CV);
                
                h_temp = (TH1D*)f_nuexsec->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_dirt", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
                h_dirt_CV = (TH1D*)h_temp->Clone();
                h_dirt_CV->Scale(_util.dirt_scale_factor / (_util.config_v.at(_util.k_Run1_Data_POT)  / POT_v.at(k_CV)) );
                h_bkg_CV->Add(h_dirt_CV);

                // HARDCODED -- NEED TO CONFIGURE A GOLBAL CALC OF FLUX AND NUM TARG
                // Scale by the flux and number of targets for the detvar CV
                h_bkg_CV->Scale(1.0/(1.40254e+10 * 4.31247e+31));
                h_bkg_CV->Scale(1.0e39);
                
                if (scale_bins)
                    h_bkg_CV->Scale(1.0, "width");
                
                h_universe.at(var).at(label).at(k_up).at(k_xsec_mcxsec_shape)->Add(h_bkg_CV, -1);
                h_universe.at(var).at(label).at(k_dn).at(k_xsec_mcxsec_shape)->Add(h_bkg_CV, -1);
                delete h_bkg_CV;
            }
        }
    
    }

    // Loop over the differential variables
    for (unsigned int var = 0; var < cv_hist_vec.size(); var++ ){
        
        
        // Loop over the types
        for (unsigned int  type = 0; type < cv_hist_vec.at(var).size(); type++ ){

            if (_util.vars.at(var) == "true_el_E" && (type != k_xsec_eff && type != k_xsec_mcxsec_smear && type != k_xsec_mcxsec_shape)) continue; // Skip the true var which doesnt make much sense

            // Get the CV histogram
            TH1D* h_CV_clone;
            
            // Use the standard CV histogram
            if (!is_detvar){
                h_CV_clone = (TH1D*)cv_hist_vec.at(var).at(type)->Clone("h_clone");
            }
            // For the beamline variations, we have a different CV
            else {
                _util.GetHist(f_nuexsec, h_CV_clone, Form( "detvar_CV/%s/h_run%s_CV_0_%s_%s", _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(type).c_str()));
                
                if (scale_bins && ((var == k_var_recoX || var == k_var_trueX) && type != k_xsec_eff ))
                        h_CV_clone->Scale(1.0, "width");
                
                h_CV_clone->SetLineWidth(2);
                h_CV_clone->SetLineColor(kBlack); 
            }

            h_CV_clone->SetMinimum(0);


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
            
            // Set the Titles
            if (type == k_xsec_mcxsec || type == k_xsec_dataxsec)
                h_CV_clone->SetTitle(_util.var_labels_xsec.at(var).c_str());
            else if (type == k_xsec_mcxsec_smear || type == k_xsec_mcxsec_shape)
                h_CV_clone->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
            else if (type == k_xsec_eff)
                h_CV_clone->SetTitle(_util.var_labels_eff.at(var).c_str());
            else
                h_CV_clone->SetTitle(_util.var_labels_events.at(var).c_str());



            // else if (type == k_xsec_eff) h_CV_clone->GetYaxis()->SetTitle("Efficiency");
            // else h_CV_clone->GetYaxis()->SetTitle("Entries");
            h_CV_clone->SetFillColorAlpha(12, 0.15);
            h_CV_clone->Draw("e2, same");
            h_CV_clone->GetXaxis()->SetTitle("");
            h_CV_clone->GetXaxis()->SetLabelSize(0);
            h_CV_clone->SetTitle(xsec_types_pretty.at(type).c_str());
            h_CV_clone->GetYaxis()->SetTitleSize(0.04);
            h_CV_clone->GetYaxis()->SetLabelSize(0.05);

            TLegend *leg;
            
            if (type != k_xsec_eff && _util.vars.at(var) != "integrated") leg = new TLegend(0.35, 0.55, 0.85, 0.85);
            else {
                
                if (is_detvar)
                    leg = new TLegend(0.35, 0.2, 0.85, 0.5);
                else
                    leg = new TLegend(0.4, 0.3, 0.9, 0.6);

            }
            leg->SetNColumns(2);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);

            leg->AddEntry(h_CV_clone,"CV (Sys.)", "lf");

            // Now draw all the universes
            for (unsigned int label = 0; label < h_universe.at(var).size(); label ++ ){
                
                // Skip the beamline variations we dont want
                if (is_detvar && (unisim_names.at(label) == "CV" || unisim_names.at(label) == "BNB_Diffusion") )
                    continue;

                SetUnisimColours(unisim_names.at(label), h_universe.at(var).at(label).at(k_up).at(type), h_universe.at(var).at(label).at(k_dn).at(type));
                

                if (single_var.at(label)){
                    
                    if(is_detvar)
                        leg->AddEntry(h_universe.at(var).at(label).at(k_up).at(type),Form("%s", var_string_pretty.at(label).c_str()), "l");
                    else
                        leg->AddEntry(h_universe.at(var).at(label).at(k_up).at(type),Form("%s", unisim_names.at(label).c_str()), "l");
                    
                    h_universe.at(var).at(label).at(k_up).at(type)->Draw("hist,same");
                }
                else {
                    if(is_detvar){
                        leg->AddEntry(h_universe.at(var).at(label).at(k_up).at(type),Form("%s +1#sigma", var_string_pretty.at(label).c_str()), "l");
                        leg->AddEntry(h_universe.at(var).at(label).at(k_dn).at(type),Form("%s -1#sigma", var_string_pretty.at(label).c_str()), "l");
                    }
                    else {
                        leg->AddEntry(h_universe.at(var).at(label).at(k_up).at(type),Form("%s +1#sigma", unisim_names.at(label).c_str()), "l");
                        leg->AddEntry(h_universe.at(var).at(label).at(k_dn).at(type),Form("%s -1#sigma", unisim_names.at(label).c_str()), "l");
                    }
                    h_universe.at(var).at(label).at(k_up).at(type)->Draw("hist,same");
                    h_universe.at(var).at(label).at(k_dn).at(type)->Draw("hist,same");
                }
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

             // Set the Titles
            if (type == k_xsec_mcxsec || type == k_xsec_dataxsec)
                h_err->SetTitle(_util.var_labels_xsec.at(var).c_str());
            else if (type == k_xsec_mcxsec_smear || type == k_xsec_mcxsec_shape)
                h_err->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
            else if (type == k_xsec_eff)
                h_err->SetTitle(_util.var_labels_eff.at(var).c_str());
            else
                h_err->SetTitle(_util.var_labels_events.at(var).c_str());

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

            // Choose what histograms to save
            if ( (var == k_var_recoX && (type == k_xsec_mcxsec || type == k_xsec_dataxsec || type == k_xsec_bkg || type == k_xsec_sig)) || 
                (var == k_var_trueX && (type == k_xsec_eff || type == k_xsec_mcxsec_smear || type == k_xsec_mcxsec_shape)) ||
                (var == k_var_integrated && (type == k_xsec_mcxsec || type == k_xsec_dataxsec || type == k_xsec_bkg || type == k_xsec_sig || type == k_xsec_eff)) 
                ){
                c->Print(Form("plots/run%s/Systematics/%s/run%s_%s_%s_%s.pdf", _util.run_period, unisim_type.c_str(), _util.run_period, unisim_type.c_str(), _util.vars.at(var).c_str(), xsec_types.at(type).c_str()));
            }

            delete c;
            delete h_CV_clone;
        }
    
    
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetUnisimColours(std::string label, TH1D* h_up, TH1D* h_dn){

    if (label == "Horn1_x" || label == "RPA"  || label == "LYRayleigh"){
        h_up->SetLineColor(95);
        h_dn->SetLineColor(95);
    }
    else if (label == "Horn_curr" || label == "CCMEC"  || label == "LYAttenuation"){
        h_up->SetLineColor(38);
        h_dn->SetLineColor(38);
    }
    else if (label == "Horn1_y" || label == "AxFFCCQE"  || label == "SCE"){
        h_up->SetLineColor(28);
        h_dn->SetLineColor(28);
    }
    else if (label == "Beam_spot" || label == "VecFFCCQE"  || label == "Recomb2"){
        h_up->SetLineColor(4);
        h_dn->SetLineColor(4);
    }
    else if (label == "Horn2_x" || label == "DecayAngMEC"  || label == "WireModX"){
        h_up->SetLineColor(kPink+1);
        h_dn->SetLineColor(kPink+1);
    }
    else if (label == "Horn2_y" || label == "ThetaDelta2Npi"  || label == "WireModYZ"){
        h_up->SetLineColor(32);
        h_dn->SetLineColor(32);
    }
    else if (label == "Horn_Water" || label == "ThetaDelta2NRad"  || label == "WireModThetaXZ"){
        h_up->SetLineColor(46);
        h_dn->SetLineColor(46);
    }
    else if (label == "Beam_shift_x" || label == "xsr_scc_Fa3"  || label == "WireModThetaYZ_withSigmaSplines"){
        h_up->SetLineColor(12);
        h_dn->SetLineColor(12);
    }
    else if (label == "Beam_shift_y" || label == "NormCCCOH"  || label == "WireModThetaYZ_withoutSigmaSplines"){
        h_up->SetLineColor(kViolet - 7);
        h_dn->SetLineColor(kViolet - 7);
    }
    else if (label == "Target_z" || label == "NormNCCOH"  || label == "WireModdEdX"){
        h_up->SetLineColor(42);
        h_dn->SetLineColor(42);
    }
    else if (label == "Decay_pipe_Bfield" || label == "LYDown" || label == "xsr_scc_Fv3"){
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
    
    std::cout << "Running code to get the reweightable variations by cut" << std::endl;
    gStyle->SetOptStat(0);

    // Load in the input file
    // Should we add more protection to this command??
    f_nuexsec = TFile::Open( Form("files/crosssec_run%s.root", _util.run_period ), "READ");

    // Define the uncertainties
    std::vector<std::tuple<std::string, int, std::string>> tuple_label = {
        std::make_tuple("RPA",              2,   "unisim"),
        std::make_tuple("CCMEC",            2,   "unisim"),
        std::make_tuple("AxFFCCQE",         2,   "unisim"),
        std::make_tuple("VecFFCCQE",        2,   "unisim"),
        std::make_tuple("DecayAngMEC",      2,   "unisim"),
        std::make_tuple("ThetaDelta2Npi",   2,   "unisim"),
        std::make_tuple("ThetaDelta2NRad",  2,   "unisim"),
        std::make_tuple("RPA_CCQE_Reduced", 2,   "unisim"),
        std::make_tuple("NormCCCOH",        2,   "unisim"),
        std::make_tuple("NormNCCOH",        2,   "unisim"),
        std::make_tuple("weightsGenie",     500, "multisim"),
        std::make_tuple("weightsReint",     1000,"multisim"),
        std::make_tuple("weightsFlux",      500, "multisim")
    };


    // Resize the vector to store the systematics
    h_cut_err.resize(tuple_label.size()); // label -- cut -- variable

    for (unsigned int label = 0; label < tuple_label.size(); label++){
        h_cut_err.at(label).resize(_util.k_cuts_MAX);
    }

    for (unsigned int label = 0; label < tuple_label.size(); label++){
        
        for (unsigned int cut = 0; cut < h_cut_err.at(label).size(); cut++){
            h_cut_err.at(label).at(cut).resize(_util.k_cut_vars_max);
        }
    }


    // Loop over cuts
    for (int cut = 0; cut < _util.k_cuts_MAX ; cut++){
        
        // Loop over the variables
        for (int var = 0; var < _util.k_cut_vars_max; var++){
            
            for (unsigned int label = 0; label < tuple_label.size(); label++){
                TH1D* h_err;
                GetCutSysUncertainty(_util.vec_hist_name.at(var), cut, std::get<0>(tuple_label.at(label)), std::get<1>(tuple_label.at(label)), std::get<2>(tuple_label.at(label)), h_err);
                
                // Save the hist in the master vector
                h_cut_err.at(label).at(cut).at(var) = (TH1D*)h_err->Clone();
                // delete h_err;
            }
        }
    }

    // Save the histograms to file
    SaveCutHistograms(tuple_label);
        
}
// -----------------------------------------------------------------------------
void SystematicsHelper::GetCutSysUncertainty(std::string histname, int cut_index, std::string label, int num_uni, std::string var_type, TH1D* &h_err){

    // if (cut_index == _util.k_cuts_MAX-1)
        // std::cout << "Saving uncertainties for: " << label  << "  " << histname << std::endl;

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
    h_err = (TH1D*)h_universe.at(k_up)->Clone(); // clone it to get the binning right

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
    for (int bin = 1; bin < h_err->GetNbinsX()+1; bin++){
        double err = std::sqrt(h_err->GetBinContent(bin)/num_uni) / h_cv->GetBinContent(bin);
        if (h_cv->GetBinContent(bin) == 0) 
            h_err->SetBinContent(bin, 0.0);
        else 
            h_err->SetBinContent(bin, err);        
    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::SaveCutHistograms(std::vector<std::tuple<std::string, int, std::string>> tuple_label){

    std::cout <<  "Now writing histograms to file..." << std::endl; 

    // ---- save the histograms into different directories inside the root file
    TFile *file_sys_var = TFile::Open(Form("files/run%s_sys_var.root", _util.run_period),"UPDATE");
    file_sys_var->cd();

    // Create subdirectory for each reweighter
    TDirectory *dir_labels[tuple_label.size()];

    // Create subdirectory for each variable
    TDirectory *dir_labels_cut[_util.k_cuts_MAX];

    // Loop over cuts
    for (int cut = 0; cut < _util.k_cuts_MAX ; cut++){

        // See if the directory already exists
        bool bool_dir = _util.GetDirectory(file_sys_var, dir_labels_cut[cut], _util.cut_dirs.at(cut).c_str());

        // If it doesnt exist then create it
        if (!bool_dir) file_sys_var->mkdir(_util.cut_dirs.at(cut).c_str());

        _util.GetDirectory(file_sys_var, dir_labels_cut[cut], _util.cut_dirs.at(cut).c_str());

        dir_labels_cut[cut]->cd();

        for (unsigned int label = 0; label < tuple_label.size(); label++){
            
            // See if the directory already exists
            bool bool_dir = _util.GetDirectory(file_sys_var, dir_labels[label], Form("%s/%s", _util.cut_dirs.at(cut).c_str(), std::get<0>(tuple_label.at(label)).c_str()));
    
            // If it doesnt exist then create it
            if (!bool_dir) file_sys_var->mkdir(Form("%s/%s", _util.cut_dirs.at(cut).c_str(), std::get<0>(tuple_label.at(label)).c_str()));

            _util.GetDirectory(file_sys_var, dir_labels[label], Form("%s/%s", _util.cut_dirs.at(cut).c_str(), std::get<0>(tuple_label.at(label)).c_str()));

            // Go into the directory
            dir_labels[label]->cd();
            
            // Loop over the variables
            for (int var = 0; var < _util.k_cut_vars_max; var++){
                h_cut_err.at(label).at(cut).at(var)->SetOption("hist");
                h_cut_err.at(label).at(cut).at(var)->Write(_util.vec_hist_name.at(var).c_str(), TObject::kOverwrite);  // write the histogram to the file
            }
            
            file_sys_var->cd();

        }

        file_sys_var->cd();
    }

    file_sys_var->cd();
    file_sys_var->Close();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SaveCutHistogramsDetVar(){

    std::cout <<  "Now writing histograms to file..." << std::endl; 

    // ---- save the histograms into different directories inside the root file
    TFile *file_sys_var = TFile::Open(Form("files/run%s_sys_var.root", _util.run_period),"UPDATE");
    file_sys_var->cd();

     // Create subdirectory for each reweighter
    TDirectory *dir_labels[k_vars_MAX+1];

    // Create subdirectory for each variable
    TDirectory *dir_labels_cut[_util.k_cuts_MAX];

    // Loop over cuts
    for (int cut = 0; cut < _util.k_cuts_MAX ; cut++){

        // See if the directory already exists
        bool bool_dir = _util.GetDirectory(file_sys_var, dir_labels_cut[cut], _util.cut_dirs.at(cut).c_str());

        // If it doesnt exist then create it
        if (!bool_dir) file_sys_var->mkdir(_util.cut_dirs.at(cut).c_str());

        _util.GetDirectory(file_sys_var, dir_labels_cut[cut], _util.cut_dirs.at(cut).c_str());

        dir_labels_cut[cut]->cd();

        for (unsigned int label = 0; label < k_vars_MAX+1; label++){
            
            // See if the directory already exists
            bool bool_dir;
            
            if (label !=  k_vars_MAX){
                bool_dir = _util.GetDirectory(file_sys_var, dir_labels[label], Form("%s/%s", _util.cut_dirs.at(cut).c_str(), var_string.at(label).c_str()));

                // If it doesnt exist then create it
                if (!bool_dir)
                    file_sys_var->mkdir(Form("%s/%s", _util.cut_dirs.at(cut).c_str(), var_string.at(label).c_str()));

                _util.GetDirectory(file_sys_var, dir_labels[label], Form("%s/%s", _util.cut_dirs.at(cut).c_str(), var_string.at(label).c_str()));
            }
            else {
                bool_dir = _util.GetDirectory(file_sys_var, dir_labels[label], Form("%s/TotalDetectorSys", _util.cut_dirs.at(cut).c_str()));

                // If it doesnt exist then create it
                if (!bool_dir)
                    file_sys_var->mkdir(Form("%s/TotalDetectorSys", _util.cut_dirs.at(cut).c_str()));

                _util.GetDirectory(file_sys_var, dir_labels[label], Form("%s/TotalDetectorSys", _util.cut_dirs.at(cut).c_str()));
            }

            // Go into the directory
            dir_labels[label]->cd();
            
            // Loop over the variables
            for (int var = 0; var < _util.k_cut_vars_max; var++){
                h_cut_err.at(label).at(cut).at(var)->SetOption("hist");
                h_cut_err.at(label).at(cut).at(var)->Write(_util.vec_hist_name.at(var).c_str(), TObject::kOverwrite);  // write the histogram to the file
            }
            
            file_sys_var->cd();

        }

        file_sys_var->cd();
    }

    file_sys_var->cd();
    file_sys_var->Close();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::MakeTotUncertaintyPlot(bool AddStatErr){


    std::vector<TH1D*> h_uncertainty;
    h_uncertainty.resize(k_ERR_MAX);

    // Loop over the variables, int, reco, true.
    for (unsigned int var = 0; var < _util.vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++){
            
            // Only look at the true efficeincy or select variables in reco space
            if ( (var == k_var_trueX && (xsec_types.at(type) == "eff" || xsec_types.at(type) == "gen_smear" || xsec_types.at(type) == "mc_xsec_smear" || xsec_types.at(type) == "mc_xsec_shape") ) || 
                 (var == k_var_recoX && (xsec_types.at(type) == "mc_xsec" || xsec_types.at(type) == "data_xsec" ||
                                             xsec_types.at(type) == "eff" || xsec_types.at(type) == "sig" || 
                                             xsec_types.at(type) == "bkg" || xsec_types.at(type) == "sel" ))){

                // Resize the vector for new loop
                if (h_uncertainty.size() == 0) 
                    h_uncertainty.resize(k_ERR_MAX);

                for (int err = 0; err < k_ERR_MAX; err++){
                    h_uncertainty.at(err) = (TH1D*)cv_hist_vec.at(var).at(type)->Clone(Form("h_clone_%s", systematic_names.at(err).c_str() ));
                }
    
                // Loop over the bins
                for (unsigned int bin = 0; bin < v_err.at(k_err_genie_uni).at(var).at(type).size(); bin++ ){

                    for (int err = 0; err < k_ERR_MAX; err++){
                        // h_uncertainty.at(err)->SetBinContent(bin+1, std::sqrt(v_err.at(err).at(var).at(type).at(bin)) );
                        
                        // In the case of MC stat error include the MC stat + response stat err in the same line
                        if (var == k_var_recoX && err == k_err_mcstats && (type == k_xsec_dataxsec || type == k_xsec_mcxsec || type == k_xsec_bkg || type == k_xsec_sig)){
                            if (type == k_xsec_dataxsec || type == k_xsec_mcxsec){
                                h_uncertainty.at(err)->SetBinContent(bin+1, 100 * std::sqrt(h_cov_v.at(var).at(type).at(err)->GetBinContent(bin+1, bin+1) + h_cov_v.at(var).at(k_xsec_mcxsec).at(k_err_stat)->GetBinContent(bin+1, bin+1)) / cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(bin+1));
                            }
                            else {
                                h_uncertainty.at(err)->SetBinContent(bin+1, 100 * std::sqrt(h_cov_v.at(var).at(type).at(err)->GetBinContent(bin+1, bin+1) + h_cov_v.at(var).at(type).at(k_err_stat)->GetBinContent(bin+1, bin+1)) / cv_hist_vec.at(var).at(type)->GetBinContent(bin+1));
                            }
                        }
                        // Add the MC stats to the MC Stat systematic to get the total stat error
                        else if (var == k_var_trueX && err == k_err_mcstats && type == k_xsec_mcxsec_shape){
                            h_uncertainty.at(err)->SetBinContent(bin+1, 100 * std::sqrt(h_cov_v.at(var).at(type).at(err)->GetBinContent(bin+1, bin+1) + h_cov_v.at(var).at(k_xsec_mcxsec_shape).at(k_err_stat)->GetBinContent(bin+1, bin+1)) / cv_hist_vec.at(var).at(k_xsec_mcxsec_shape)->GetBinContent(bin+1));
                        }
                        // Plot data stat error in both cases
                        else if (var == k_var_recoX && err == k_err_stat && (type == k_xsec_dataxsec || type == k_xsec_mcxsec)){
                            h_uncertainty.at(err)->SetBinContent(bin+1, 100 * std::sqrt(h_cov_v.at(var).at(k_xsec_dataxsec).at(k_err_stat)->GetBinContent(bin+1, bin+1)) / cv_hist_vec.at(var).at(k_xsec_dataxsec)->GetBinContent(bin+1));
                        }
                        // Plot data stat
                        else if (var == k_var_trueX && err == k_err_stat && (type == k_xsec_mcxsec_shape)){
                            h_uncertainty.at(err)->SetBinContent(bin+1, 100 * std::sqrt(h_cov_v.at(var).at(k_xsec_dataxsec).at(k_err_stat)->GetBinContent(bin+1, bin+1)) / cv_hist_vec.at(var).at(k_xsec_dataxsec)->GetBinContent(bin+1));
                        }
                        else {
                            h_uncertainty.at(err)->SetBinContent(bin+1, 100 * std::sqrt(h_cov_v.at(var).at(type).at(err)->GetBinContent(bin+1, bin+1)) / cv_hist_vec.at(var).at(type)->GetBinContent(bin+1));
                        }

                        h_uncertainty.at(err)->SetLineWidth(2);
                    }
                    
                }

                // Choose if we want to draw only the sys error or sys + stat error
                if (AddStatErr && (type == k_xsec_dataxsec || type == k_xsec_mcxsec || type == k_xsec_mcxsec_shape)){
                    h_uncertainty.at(k_err_tot)->SetLineColor(kBlack);
                    h_uncertainty.at(k_err_stat)->SetLineColor(kBlack);
                    h_uncertainty.at(k_err_stat)->SetLineStyle(2);
                }
                else { 
                    h_uncertainty.at(k_err_sys)->SetLineColor(kBlack);
                }

                
                h_uncertainty.at(k_err_genie_uni)->SetLineColor(95);
                h_uncertainty.at(k_err_genie_multi)->SetLineColor(kPink+1);
                h_uncertainty.at(k_err_beamline)->SetLineColor(28);
                h_uncertainty.at(k_err_hp)->SetLineColor(4);
                h_uncertainty.at(k_err_reint)->SetLineColor(kViolet-1);
                h_uncertainty.at(k_err_detvar)->SetLineColor(32);
                h_uncertainty.at(k_err_dirt)->SetLineColor(46);
                h_uncertainty.at(k_err_pot)->SetLineColor(42);
                h_uncertainty.at(k_err_ext)->SetLineColor(kSpring-1);
                h_uncertainty.at(k_err_mcstats)->SetLineColor(7);

                TCanvas *c = new TCanvas("c", "c", 500, 500);
                c->SetLeftMargin(0.15);
                c->SetBottomMargin(0.13);

                if (type == k_xsec_mcxsec_shape){
                    h_uncertainty.at(k_err_tot)->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
                    h_uncertainty.at(k_err_sys)->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
                }
                
                if (AddStatErr && (type == k_xsec_dataxsec || type == k_xsec_mcxsec || type == k_xsec_mcxsec_shape)){
                    h_uncertainty.at(k_err_tot)->SetTitle(xsec_types_pretty.at(type).c_str());
                    h_uncertainty.at(k_err_tot)->GetYaxis()->SetTitle("Uncertainty [%]");
                    
                    if (std::string(_util.xsec_var) == "elec_E")
                        h_uncertainty.at(k_err_tot)->GetYaxis()->SetRangeUser(0, 250);
                    else 
                        h_uncertainty.at(k_err_tot)->GetYaxis()->SetRangeUser(0, 100);

                    // if (type == k_xsec_mcxsec_shape)
                    //     h_uncertainty.at(k_err_tot)->GetYaxis()->SetRangeUser(0, 250);
                    
                    h_uncertainty.at(k_err_tot)->Draw("hist,same, text00");
                    h_uncertainty.at(k_err_stat)->Draw("hist,same");
                }
                else {
                    h_uncertainty.at(k_err_sys)->SetTitle(xsec_types_pretty.at(type).c_str());
                    h_uncertainty.at(k_err_sys)->GetYaxis()->SetTitle("Uncertainty [%]");
                    
                    if (std::string(_util.xsec_var) == "elec_E")
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 150);
                    else 
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 100);
                    
                    h_uncertainty.at(k_err_sys)->Draw("hist,same, text00");
                }

                if (type == k_xsec_bkg){
                    if (std::string(_util.xsec_var) == "elec_E")
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 240);
                    else
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 120);
                }

                if (type == k_xsec_mcxsec_smear){
                    if (std::string(_util.xsec_var) == "elec_E")
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 100);
                    else
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 50);
                }

                if (type == k_xsec_sig){
                    if (std::string(_util.xsec_var) == "elec_E")
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 90);
                    else
                        h_uncertainty.at(k_err_sys)->GetYaxis()->SetRangeUser(0, 60);
                }
                
                h_uncertainty.at(k_err_genie_uni)->Draw("hist,same");
                h_uncertainty.at(k_err_genie_multi)->Draw("hist,same");
                h_uncertainty.at(k_err_beamline)->Draw("hist,same");
                h_uncertainty.at(k_err_hp)->Draw("hist,same");
                h_uncertainty.at(k_err_reint)->Draw("hist,same");
                h_uncertainty.at(k_err_detvar)->Draw("hist,same");
                h_uncertainty.at(k_err_dirt)->Draw("hist,same");
                h_uncertainty.at(k_err_pot)->Draw("hist,same");
                h_uncertainty.at(k_err_ext)->Draw("hist,same");
                h_uncertainty.at(k_err_mcstats)->Draw("hist,same");

                TLegend *leg = new TLegend(0.18, 0.55, 0.88, 0.85);
                leg->SetNColumns(2);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                if (AddStatErr && (type == k_xsec_dataxsec || type == k_xsec_mcxsec || type == k_xsec_mcxsec_shape)){
                    leg->AddEntry(h_uncertainty.at(k_err_tot),        "Total Stat. + Sys.", "l");
                    leg->AddEntry(h_uncertainty.at(k_err_stat),       "Data Stat.", "l");
                }
                else {
                    leg->AddEntry(h_uncertainty.at(k_err_sys),        "Total Sys.", "l");
                }
                leg->AddEntry(h_uncertainty.at(k_err_hp),         "Hadron Production", "l");
                leg->AddEntry(h_uncertainty.at(k_err_detvar),     "Detector", "l");
                leg->AddEntry(h_uncertainty.at(k_err_genie_multi),"GENIE Multisim", "l");
                leg->AddEntry(h_uncertainty.at(k_err_genie_uni),  "GENIE Unisim", "l");
                leg->AddEntry(h_uncertainty.at(k_err_beamline),   "Beamline Geometry", "l");
                leg->AddEntry(h_uncertainty.at(k_err_reint),      "Geant4 Reinteractions", "l");
                leg->AddEntry(h_uncertainty.at(k_err_pot),        "POT Counting", "l");
                leg->AddEntry(h_uncertainty.at(k_err_dirt),       "Dirt", "l");
                leg->AddEntry(h_uncertainty.at(k_err_ext),        "EXT Norm", "l");
                leg->AddEntry(h_uncertainty.at(k_err_mcstats),    "MC Stat", "l");

                leg->Draw();

                if (AddStatErr && (type == k_xsec_dataxsec || type == k_xsec_mcxsec || type == k_xsec_mcxsec_shape)){
                    c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_%s_%s_tot_uncertainty.pdf", _util.run_period, _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(type).c_str()));
                }
                else {
                    c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_%s_%s_tot_sys_uncertainty.pdf", _util.run_period, _util.vars.at(var).c_str(), _util.run_period, _util.vars.at(var).c_str(), xsec_types.at(type).c_str()));
                }

                h_uncertainty.clear();
                delete c;

            }

        }
        
        
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialseCovarianceVector(){

    // Initialise the covariance matrices -- this needs to be vectorized
    int n_bins = cv_hist_vec.at(k_var_recoX).at(0)->GetNbinsX();
    
    // Resize to the number of variables
    h_cov_v.resize(_util.vars.size());

    // Loop over _util.vars and resize to the number of types
    for (unsigned int var = 0; var < h_cov_v.size(); var++) {
        h_cov_v.at(var).resize(xsec_types.size());
    }
    
    // Loop over _util.vars
    for (unsigned int var = 0; var < h_cov_v.size(); var++) {
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++) {
            // Resize to the number of systematics 
            h_cov_v.at(var).at(type).resize(k_ERR_MAX);
        }
    
    }

    // Do the looping again and create the histograms
    for (unsigned int var = 0; var < h_cov_v.size(); var++) {
        for (unsigned int type = 0; type < h_cov_v.at(var).size(); type++) {
            for (unsigned int cov = 0; cov < h_cov_v.at(var).at(type).size(); cov++){
                
                if (var == k_var_integrated){
                    h_cov_v.at(var).at(type).at(cov) = new TH2D(Form("h_cov_%s_%s_%s", _util.vars.at(var).c_str(), xsec_types.at(type).c_str(), systematic_names.at(cov).c_str()), "Covariance Matrix ;Bin i; Bin j", 1, 1, 2, 1, 1, 2);
                }
                else { 
                    h_cov_v.at(var).at(type).at(cov) = new TH2D(Form("h_cov_%s_%s_%s", _util.vars.at(var).c_str(), xsec_types.at(type).c_str(), systematic_names.at(cov).c_str()), "Covariance Matrix ;Bin i; Bin j", n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
                }
                
            }
        }
    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::ExportResult(TFile* f){

    std::cout << _util.blue << "Saving Results to File!!!" << _util.reset << std::endl;

    // First get the CV/response matrix to store to the file
    TH2D* h_response;
    TH1D* h_mcxsec_fine;

    int _var = k_var_recoX;
    std::string folder_name = "electron_energy";
    
    // MCC8 so get the smearing matrix
    if (std::string(_util.xsec_smear_mode) == "mcc8"){
        h_response = (TH2D*)f->Get(Form("CV/%s/h_run%s_CV_0_smearing",_util.vars.at(k_var_recoX).c_str(),_util.run_period));
    }
    // Other modes we need the response matrix
    else {
        h_response = (TH2D*)f->Get(Form("CV/%s/h_run%s_CV_0_smearing",_util.vars.at(k_var_trueX).c_str(),_util.run_period));
        h_mcxsec_fine = (TH1D*)f->Get(Form("CV/%s/h_run%s_CV_0_%s_mc_xsec_fine",_util.vars.at(k_var_trueX).c_str(),_util.run_period, _util.vars.at(k_var_trueX).c_str()));
        _var = k_var_trueX;
    }

    // Choose whichj folder to put the stuff in
    if (std::string(_util.xsec_var) == "elec_ang"){
        folder_name = "elec_ang";
    }
    else if (std::string(_util.xsec_var) == "elec_cang"){
        folder_name = "elec_cang";
    }
    else{
        folder_name = "elec_E";
    }

    // Now we have the histograms we need, lets open a new file to store the systematics
    TFile *f_sys_out;
    if (std::string(_util.xsec_bin_mode) == "ratio"){
        f_sys_out = TFile::Open(Form("files/xsec_result_run%s_ratio.root", _util.run_period), "UPDATE");
    }
    else {
        f_sys_out = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "UPDATE");
    }
    
    f_sys_out->cd();

    // Create subdirectory for each reweighter
    TDirectory *dir_mode;

    // Create subdirectory for each variable
    TDirectory *dir_var;

    // See if the directory already exists
    bool bool_dir = _util.GetDirectory(f_sys_out, dir_var, folder_name.c_str());

    // If it doesnt exist then create it
    if (!bool_dir) f_sys_out->mkdir(folder_name.c_str());

    _util.GetDirectory(f_sys_out, dir_var, folder_name.c_str());

    dir_var->cd();

    // See if the directory already exists
    bool_dir = _util.GetDirectory(f_sys_out, dir_mode, Form("%s/%s", folder_name.c_str(), _util.xsec_smear_mode));

    // If it doesnt exist then create it
    if (!bool_dir) f_sys_out->mkdir(Form("%s/%s", folder_name.c_str(), _util.xsec_smear_mode));

    _util.GetDirectory(f_sys_out, dir_mode, Form("%s/%s", folder_name.c_str(), _util.xsec_smear_mode));

    // Go into the directory
    dir_mode->cd();
    
    // MCC8 mode, so get the smearing matrix
    if (std::string(_util.xsec_smear_mode) == "mcc8"){

        // Smearing Matrix  ---------------------------------
        h_response->SetOption("col,text00");
        h_response->Write("h_smear", TObject::kOverwrite);

        // MC XSec Covariance Matrix  ---------------------------------
        h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_tot)->SetOption("col");
        h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_tot)->Write("h_cov_tot_mcxsec_reco", TObject::kOverwrite);

        // MC XSec Reco  ---------------------------------
        cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec)->SetOption("hist");
        cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec)->Write("h_mc_xsec_reco", TObject::kOverwrite);

        // MC Efficiency Smear  ---------------------------------
        cv_hist_vec.at(k_var_recoX).at(k_xsec_eff)->SetOption("hist");
        cv_hist_vec.at(k_var_recoX).at(k_xsec_eff)->Write("h_mc_eff_reco", TObject::kOverwrite);

        // MC Correlation Matrix  ---------------------------------
        TH2D* cor;
        int n_bins = cv_hist_vec.at(k_var_recoX).at(0)->GetNbinsX();
        cor  = new TH2D("h_cor",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the correlation matrix
        _util.CalcCorrelation(cv_hist_vec.at(k_var_recoX).at(k_xsec_mcxsec), h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_tot), cor);

        cor->SetOption("colz");
        cor->Write("h_corr_tot_mcxsec_reco", TObject::kOverwrite);

    }
    // Other modes we need the response matrix
    else if (std::string(_util.xsec_smear_mode) == "er") {
        
        // Response Matrix  ---------------------------------
        h_response->SetOption("col,text00");
        h_response->Write("h_response", TObject::kOverwrite);

        // MC XSec Covariance Matrix  ---------------------------------
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot)->SetOption("col");
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot)->Write("h_cov_tot_mcxsec_reco", TObject::kOverwrite);

        // MC XSec Sys Covariance Matrix  ---------------------------------
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_sys)->SetOption("col");
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_sys)->Write("h_cov_sys_mcxsec_reco", TObject::kOverwrite);

        // MC XSec Covariance Matrix  ---------------------------------
        TH2D* h_cov_xsec_sys_tot = (TH2D*)h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_genie_multi)->Clone();
        h_cov_xsec_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_genie_uni));
        h_cov_xsec_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_mcstats));
        // h_cov_xsec_sys_tot->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec_smear).at(k_err_genie_multi));
        // h_cov_xsec_sys_tot->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec_smear).at(k_err_genie_uni));
        // h_cov_xsec_sys_tot->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec_smear).at(k_err_mcstats));

        h_cov_xsec_sys_tot->SetOption("col");
        h_cov_xsec_sys_tot->Write("h_cov_xsec_sys_mcxsec_reco", TObject::kOverwrite);

        // MC XSec Flux Covariance Matrix  ---------------------------------
        TH2D* h_cov_flux_sys_tot = (TH2D*)h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_hp)->Clone();
        h_cov_flux_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_beamline));
        h_cov_flux_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_mcstats));
        // h_cov_flux_sys_tot->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec_smear).at(k_err_hp));
        // h_cov_flux_sys_tot->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec_smear).at(k_err_beamline));
        // h_cov_flux_sys_tot->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_stat));
        h_cov_flux_sys_tot->SetOption("col");
        h_cov_flux_sys_tot->Write("h_cov_flux_sys_mcxsec_reco", TObject::kOverwrite);

        // MC XSec Genie All Covariance Matrix  ---------------------------------
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_genie_multi)->SetOption("col");
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_genie_multi)->Write("h_cov_genie_multi_mcxsec_reco", TObject::kOverwrite);

        // MC XSec Smear MC Stats Covariance Matrix  ---------------------------------
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_mcstats)->SetOption("col");
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_mcstats)->Write("h_cov_tot_mcxsec_smear_true", TObject::kOverwrite);
    
        // MC XSec Reco  ---------------------------------
        cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape)->SetOption("hist");
        cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape)->Write("h_mc_xsec_reco", TObject::kOverwrite);

        // MC Correlation Matrix  ---------------------------------
        TH2D* cor;
        int n_bins = cv_hist_vec.at(k_var_recoX).at(0)->GetNbinsX();
        cor  = new TH2D("h_cor",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the correlation matrix
        _util.CalcCorrelation(cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape), h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot), cor);

        cor->SetOption("colz");
        cor->Write("h_corr_tot_mcxsec_reco", TObject::kOverwrite);

    }
    // Weiner Mode
    else {
        std::vector<double> temp_bins = { 0.23, 0.41, 0.65, 0.94, 1.35, 1.87, 2.32, 6.0};
        double* edges = &temp_bins[0]; // Cast to an array 

        std::vector<double> bins = { 0.0, 0.23, 0.41, 0.65, 0.94, 1.35, 1.87, 2.32, 6.0};
        double* edges2 = &bins[0]; // Cast to an array 

        std::vector<double> bins_fine = _util.true_shr_bins;
        double* edges_fine = &bins_fine[0]; // Cast to an array 

        // Response Matrix  ---------------------------------
        h_response->SetOption("col,text00");
        
        // Flip the response matrix x and y axes
        // TH2D* h_smear = new TH2D("h_response","; Leading Shower Energy [GeV]; True e#lower[-0.5]{-} + e^{+} Energy [GeV]", bins.size(), edges2, bins_fine.size(), edges_fine);
        TH2D* h_smear = (TH2D*)h_response->Clone();

        // Loop over rows
        for (int row=1; row<h_response->GetXaxis()->GetNbins()+1; row++) {

            for (int col=1; col<h_response->GetYaxis()->GetNbins()+1; col++){
                h_smear->SetBinContent(col, row, h_response->GetBinContent(row, col));          
            }
        } 
        
        if (std::string(_util.xsec_var) == "elec_ang"){
            h_smear->SetTitle("; #beta^{reco}_{e#lower[-0.5]{-} + e^{+}} [deg]; #beta^{true}_{e#lower[-0.5]{-} + e^{+}} [deg]");
        }
        else if (std::string(_util.xsec_var) == "elec_cang"){
            h_smear->SetTitle("; cos(#beta)^{reco}_{e#lower[-0.5]{-} + e^{+}} [deg]; #beta^{true}_{e#lower[-0.5]{-} + e^{+}} [deg]");
        }
        else{
            h_smear->SetTitle("; E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV]; E^{true}_{e#lower[-0.5]{-} + e^{+}} [GeV]");
        }
        
        
        h_smear->Write("h_response", TObject::kOverwrite);

        // Data XSec Covariance Matrix ---------------------------------
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot)->SetOption("colz");
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot)->Write("h_cov_tot_dataxsec_reco", TObject::kOverwrite);

        // Data XSec Stat Covariance Matrix ---------------------------------
        h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(k_err_stat)->SetOption("colz");
        h_cov_v.at(k_var_recoX).at(k_xsec_dataxsec).at(k_err_stat)->Write("h_cov_stat_dataxsec_reco", TObject::kOverwrite);

        // Data XSec Sys Covariance Matrix ---------------------------------
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_sys)->SetOption("colz");
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_sys)->Write("h_cov_sys_dataxsec_reco", TObject::kOverwrite);

        // Data XSec Covariance Matrix  ---------------------------------
        TH2D* h_cov_xsec_sys_tot = (TH2D*)h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_genie_multi)->Clone();
        h_cov_xsec_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_genie_uni));
        h_cov_xsec_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_mcstats));
        
        // // Create a temp matrix to gather all the MC smear and stat before adding to the data
        // TH2D* h_cov_xsec_sys_tot_temp = (TH2D*)h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_stat)->Clone();
        // h_cov_xsec_sys_tot_temp->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_smear).at(k_err_genie_multi));
        // h_cov_xsec_sys_tot_temp->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_smear).at(k_err_genie_uni));
        // h_cov_xsec_sys_tot_temp->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_smear).at(k_err_mcstats));

        h_cov_xsec_sys_tot->SetOption("col");
        h_cov_xsec_sys_tot->Write("h_cov_xsec_sys_dataxsec_reco", TObject::kOverwrite);

        // Data Flux Covariance Matrix  ---------------------------------
        TH2D* h_cov_flux_sys_tot = (TH2D*)h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_hp)->Clone();
        h_cov_flux_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_beamline));
        h_cov_flux_sys_tot->Add(h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_mcstats));
        
        // Create a temp matrix to gather all the MC smear and stat before adding to the data
        // TH2D* h_cov_flux_sys_tot_temp = (TH2D*)h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec).at(k_err_stat)->Clone();
        // h_cov_flux_sys_tot_temp->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec_smear).at(k_err_hp));
        // h_cov_flux_sys_tot_temp->Add(h_cov_v.at(k_var_recoX).at(k_xsec_mcxsec_smear).at(k_err_beamline));
        h_cov_flux_sys_tot->SetOption("col");
        h_cov_flux_sys_tot->Write("h_cov_flux_sys_dataxsec_reco", TObject::kOverwrite);
    
        // MC XSec True  ---------------------------------
        cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec)->SetOption("hist");
        cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec)->Write("h_mc_xsec_true", TObject::kOverwrite);

        // MC Efficiency True  ---------------------------------
        cv_hist_vec.at(k_var_trueX).at(k_xsec_eff)->SetOption("hist");
        cv_hist_vec.at(k_var_trueX).at(k_xsec_eff)->Write("h_mc_eff_true", TObject::kOverwrite);

        // MC XSec Covariance Matrix ---------------------------------
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot)->SetOption("colz");
        h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot)->Write("h_cov_tot_mcxsec_reco", TObject::kOverwrite);


        // Data Correlation Matrix  ---------------------------------
        TH2D* cor;
        int n_bins = cv_hist_vec.at(k_var_recoX).at(0)->GetNbinsX();
        cor  = new TH2D("h_cor",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the correlation matrix
        _util.CalcCorrelation(cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec_shape), h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot), cor);

        cor->SetOption("colz");
        cor->Write("h_corr_tot_dataxsec_reco", TObject::kOverwrite);

        // Unfolded result
        _wSVD.DoUnfolding(2, 0, cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec), cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec), h_smear, h_cov_v.at(k_var_trueX).at(k_xsec_mcxsec_shape).at(k_err_tot));
        _wSVD.CompareModel(cv_hist_vec.at(k_var_trueX).at(k_xsec_mcxsec));
        _wSVD.unf->Write("h_data_xsec_unfolded", TObject::kOverwrite);
        _wSVD.unfcov->Write("h_data_cov_tot_unfolded", TObject::kOverwrite);
        _wSVD.smear->Write("h_ac", TObject::kOverwrite);
    }
    

    // Data XSec Reco (Stat Only)  ---------------------------------
    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetOption("E1,X0");
    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetLineColor(kRed+2);
    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetTitle(_util.var_labels_xsec.at(k_var_recoX).c_str());
    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->Write("h_data_xsec_stat_reco", TObject::kOverwrite);

    // Data XSec Reco (Sys Only)  ---------------------------------
    for (int bin = 0; bin < cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->GetNbinsX(); bin++){
        cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetBinError(bin+1, 0.01*std::sqrt(v_err.at(k_err_sys).at(k_var_trueX).at(k_xsec_mcxsec_shape).at(bin)) * cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->GetBinContent(bin+1));
    }

    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetOption("E1,X0");
    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->Write("h_data_xsec_sys_reco", TObject::kOverwrite);


    // Data XSec Reco (Stat + Sys)  ---------------------------------
    for (int bin = 0; bin < cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->GetNbinsX(); bin++){
        cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetBinError(bin+1, 0.01*std::sqrt(v_err.at(k_err_sys).at(k_var_trueX).at(k_xsec_mcxsec_shape).at(bin) + v_err.at(k_err_stat).at(k_var_recoX).at(k_xsec_dataxsec).at(bin)) * cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->GetBinContent(bin+1));
    }

    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetOption("E1,X0");
    cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->Write("h_data_xsec_stat_sys_reco", TObject::kOverwrite);


    // Close the file ---------------------------------
    f_sys_out->cd();
    f_sys_out->Close();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::ExportTotalCrossSectionResult(){

    std::cout << _util.blue << "Saving Total Cross Section Results to File!!!" << _util.reset << std::endl;

    std::string folder_name = "total";
    
    // Now we have the histograms we need, lets open a new file to store the systematics
    TFile *f_sys_out = TFile::Open(Form("files/xsec_result_run%s.root", _util.run_period), "UPDATE");
    f_sys_out->cd();

    // Create subdirectory for each reweighter
    TDirectory *dir_mode;

    // Create subdirectory for each variable
    TDirectory *dir_var;

    // See if the directory already exists
    bool bool_dir = _util.GetDirectory(f_sys_out, dir_var, folder_name.c_str());

    // If it doesnt exist then create it
    if (!bool_dir) f_sys_out->mkdir(folder_name.c_str());

    _util.GetDirectory(f_sys_out, dir_var, folder_name.c_str());

    dir_var->cd();

    // Now save the histograms
    // Data XSec Reco (Stat Only)  ---------------------------------
    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->SetOption("E1,X0");
    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->SetLineColor(kBlack);
    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->SetTitle(_util.var_labels_xsec.at(k_var_integrated).c_str());
    // Data XSec Reco (Sys Only)  ---------------------------------
    for (int bin = 0; bin < cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->GetNbinsX(); bin++){
        cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetBinError(bin+1, 0.01*std::sqrt(v_err.at(k_err_stat).at(k_var_integrated).at(k_xsec_mcxsec).at(bin)) * cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->GetBinContent(bin+1));
    }

    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->Write("h_data_xsec_stat_reco", TObject::kOverwrite);

    // Data XSec Reco (Sys Only)  ---------------------------------
    for (int bin = 0; bin < cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->GetNbinsX(); bin++){
        cv_hist_vec.at(k_var_recoX).at(k_xsec_dataxsec)->SetBinError(bin+1, 0.01*std::sqrt(v_err.at(k_err_sys).at(k_var_integrated).at(k_xsec_mcxsec).at(bin)) * cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->GetBinContent(bin+1));
    }

    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->SetOption("E1,X0");
    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->Write("h_data_xsec_sys_reco", TObject::kOverwrite);


    // Data XSec Reco (Stat + Sys)  ---------------------------------
    for (int bin = 0; bin < cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->GetNbinsX(); bin++){
        cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->SetBinError(bin+1, 0.01*std::sqrt(v_err.at(k_err_sys).at(k_var_integrated).at(k_xsec_mcxsec).at(bin) + v_err.at(k_err_stat).at(k_var_integrated).at(k_xsec_dataxsec).at(bin)) * cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->GetBinContent(bin+1));
    }

    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->SetOption("E1,X0");
    cv_hist_vec.at(k_var_integrated).at(k_xsec_dataxsec)->Write("h_data_xsec_stat_sys_reco", TObject::kOverwrite);

    // Close the file ---------------------------------
    f_sys_out->cd();
    f_sys_out->Close();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::MakedEdxPaperPlot(){

    // ------------------------------------------------------------
    // Some initial configurations and work around fixes

    int cut = _util.k_unselected;
    const char* x_axis_name = "Leading Shower dE/dx [MeV/cm]";
    bool plotdata = true;
    
    
    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_diff; // The vector of histogram differentes between CV and the vatiation (variation-CV)
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    hist.resize(k_vars_MAX);
    hist_diff.resize(k_vars_MAX);
    hist_ratio.resize(k_vars_MAX);

    // Get the data histogram
    TFile *f_data;
    TH1D *h_data, *h_dirt, *h_ext;

    if (plotdata){
        f_data= TFile::Open("files/nuexsec_run1_merged.root", "READ");
        _util.GetHist(f_data, h_data, Form("Stack/%s/%s/%s_%s_%s", _util.cut_dirs.at(cut).c_str(), "data", "h_reco_shr_tkfit_dedx_max_tune", _util.cut_dirs.at(cut).c_str(), "data"));
        _util.GetHist(f_data, h_dirt, Form("Stack/%s/%s/%s_%s_%s", _util.cut_dirs.at(cut).c_str(), "dirt", "h_reco_shr_tkfit_dedx_max_tune", _util.cut_dirs.at(cut).c_str(), "dirt"));
        _util.GetHist(f_data, h_ext, Form("Stack/%s/%s/%s_%s_%s",  _util.cut_dirs.at(cut).c_str(), "ext",  "h_reco_shr_tkfit_dedx_max_tune", _util.cut_dirs.at(cut).c_str(), "ext"));
        h_data->SetDirectory(0);
        h_dirt->SetDirectory(0);
        h_ext->SetDirectory(0);
        f_data->Close();
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.5);
        h_data->SetLineColor(kBlack);
        h_dirt->Scale(_util.dirt_scale_factor);
        h_ext->Scale(_util.ext_scale_factor);
    }

    TCanvas * c      = new TCanvas("","",500,500);
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    _util.SetTPadOptions(topPad, bottomPad );

    // ------------------------------------------------------------
    // Loop over the variations and get the histograms
    
    for (unsigned int k=0; k < f_vars.size(); k++){
        
        // Loop over the classifications and get the histograms
        for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

            // Only want the MC piece -- may want to add in dirt too? -- will need to separately scale that hitogram though
            if ( i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt ) continue;

            // Get all the MC histograms and add them together
            TH1D *h_temp;

            _util.GetHist(f_vars.at(k), h_temp, Form("Stack/%s/%s/%s_%s_%s", _util.cut_dirs.at(cut).c_str(), _util.classification_dirs.at(i).c_str(),  "h_reco_shr_tkfit_dedx_max_tune", _util.cut_dirs.at(cut).c_str(), _util.classification_dirs.at(i).c_str()));

            // First case so clone the histogram
            if (i == 0) hist.at(k) = (TH1D*) h_temp->Clone("h_sum_hist");
            else hist.at(k)->Add(h_temp, 1);
        }
        
    }

    // Legend
    // on top of the topCanvs to avoid overlaping the plot
    TLegend *leg = new TLegend(0.1686747,0.7233083,0.8795181,0.8406015,NULL,"brNDC");
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // ------------------------------------------------------------
    // Now scale the histograms to POT
    // and draw the histograms on the top pad

    // variable used to track the max bin content to later scale the plot to avoid cutting information
    Double_t max_bin = 0.; 

    // Clone a histogram to plot the CV error as a grey band
    TH1D* h_error_hist_data = (TH1D*)h_data->Clone();
    h_error_hist = (TH1D*) hist.at(k_CV)->Clone("h_error_hist");
    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->SetLineWidth(2);
    h_error_hist->SetLineColor(kBlack);
    h_error_hist->Scale(_util.mc_scale_factor*3.0754529);
    h_error_hist->Add(h_dirt, 1.0);
    h_error_hist->Add(h_ext, 1.0);

    // Genie Unisim
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "RPA", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "CCMEC", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "AxFFCCQE", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "VecFFCCQE", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "DecayAngMEC", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "ThetaDelta2Npi", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "ThetaDelta2NRad", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "NormCCCOH", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "NormNCCOH", "MC");
    
    // Beamline
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Horn1_x", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Horn_curr", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Horn1_y", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Beam_spot", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Horn2_x", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Horn2_y", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Horn_Water", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Beam_shift_x", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Beam_shift_y", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Target_z", "MC");
            
    // Other
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "weightsFlux", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "weightsGenie", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "weightsReint", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "POT",  "Stack");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "Dirt", "Dirt");

    // Clone a histogram to plot the CV error as a grey band
    TH1D* h_error_hist_noDetvar = (TH1D*) h_error_hist->Clone("h_error_hist_nodetvar");
    h_error_hist_noDetvar->SetFillColorAlpha(kRed+2, 0.15);

    // Individual detector systematics
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "WireModX", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "WireModYZ", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "WireModThetaXZ", "MC");
    AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max_tune", "Moliere_Avg", "WireModThetaYZ_withoutSigmaSplines", "MC");

    // AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max", "Moliere_Avg", "LYRayleigh", "MC");
    // AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max", "Moliere_Avg", "SCE", "MC");
    // AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max", "Moliere_Avg", "Recomb2", "MC");
    // AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max", "Moliere_Avg", "WireModThetaYZ_withSigmaSplines", "MC");
    // AddSysUncertainty(h_error_hist, h_ext, h_dirt, "h_reco_shr_tkfit_dedx_max", "Moliere_Avg", "WireModdEdX", "MC");




    for (unsigned int y=0; y < hist.size(); y++ ){
        double scale_fact = POT_v.at(k_CV) / POT_v.at(y);
        // std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        if (plotdata){
            hist.at(y)->Scale(_util.mc_scale_factor*3.0754529); // extra factor to scale to main mc pot rather than detvar -- hack
            hist.at(y)->Add(h_dirt, 1.0);
            hist.at(y)->Add(h_ext, 1.0);
        }

        // Save clones of the histograms for doing the ratios
        hist_ratio.at(y) = (TH1D*) hist.at(y)->Clone(Form("h_ratio_%s", var_string.at(y).c_str()));
        hist_ratio.at(y)->Divide(hist.at(k_CV));

        // Set the customisation of the histogram
        SetVariationProperties(hist.at(y), y);
        SetVariationProperties(hist_ratio.at(y), y);

        // Draw the histograms
        if (y == k_CV) {
            leg->AddEntry(h_error_hist, "CV (Tot Unc.)", "lf");
            
            if (hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()) > max_bin)
                max_bin = hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()); // for scale purposes
        
        }
        else {
            hist.at(y)->SetLineColor(var_string_pretty_color.at(y));  // change color of the histogram
            hist_ratio.at(y)->SetLineColor(var_string_pretty_color.at(y));
            leg->AddEntry(hist.at(y), var_string_pretty.at(y).c_str(), "l"); // add histogram to legend
            
            if (hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()) > max_bin)
                max_bin = hist.at(y)->GetBinContent(hist.at(y)->GetMaximumBin()); // for scale purposes
        }
    }
    

    if (plotdata)
         leg->AddEntry(h_data, "Beam-On", "lep"); // add histogram to legend

    // -----------------------------------------------------------------
    // Drawing histograms on top pad

    // setting hist config to the first one that we drawn
    hist.at(0)->GetYaxis()->SetRangeUser(0,max_bin*1.4);
    hist.at(0)->GetYaxis()->SetTitle("Entries");
    hist.at(0)->GetXaxis()->SetLabelSize(0);
    hist.at(0)->GetYaxis()->SetTitleSize(0.05);
    hist.at(0)->GetYaxis()->SetLabelSize(0.05);

    // drawing histograms
    for (unsigned int y=0; y < hist.size(); y++ ) {
    
            hist.at(0)->GetXaxis()->SetRangeUser(0,7);
            hist.at(y)->Draw("hist, same");
            //if (y == k_CV) h_error_hist->Draw("E2, same");
    }
    
    // -----------------------------------------------------------------
    
    h_error_hist->Draw("E2, same");

    TH1D * h_error_hist_ratio = (TH1D*)h_error_hist->Clone();
    TH1D * h_error_hist_noDetvar_ratio = (TH1D*)h_error_hist->Clone();
    // Take the ratio of the error histograms
    for(int bin = 1; bin <= h_error_hist->GetNbinsX(); bin++){
        h_error_hist_ratio->SetBinContent( bin , 1.0);
        h_error_hist_ratio->SetBinError( bin , h_error_hist->GetBinError(bin)/h_error_hist->GetBinContent(bin) );
        h_error_hist_noDetvar_ratio->SetBinContent( bin , 1.0);
        h_error_hist_noDetvar_ratio->SetBinError( bin , h_error_hist_noDetvar->GetBinError(bin)/h_error_hist->GetBinContent(bin) );
        
        h_error_hist_data->SetBinContent( bin , h_error_hist_data->GetBinContent(bin)/h_error_hist->GetBinContent(bin) );
        h_error_hist_data->SetBinError( bin , h_error_hist_data->GetBinError(bin)/h_error_hist->GetBinContent(bin) );
    }

    // h_error_hist_data->Divide(h_error_hist);

    // h_error_hist_ratio = (TH1D*) h_error_hist->Clone("h_error_hist_rat");
    // h_error_hist_ratio->Divide();

    // drawing CV again to make sure it is on top of everything else
    hist.at(k_CV)->Draw("hist, same");

    if (plotdata){
        h_data->Draw("same PE");
    }

    // drawing legend
    leg->Draw();


    // -----------------------------------------------------------------
    // Drawing the total detector systematic uncertainty on the bottom pad

    bottomPad->cd();

    hist_ratio.at(0)->SetLineWidth(2);
    hist_ratio.at(0)->GetXaxis()->SetLabelSize(15); // 12
    hist_ratio.at(0)->GetXaxis()->SetLabelFont(43); 
    hist_ratio.at(0)->GetYaxis()->SetLabelSize(11);
    hist_ratio.at(0)->GetYaxis()->SetLabelFont(43);
    hist_ratio.at(0)->GetXaxis()->SetTitleOffset(3.2); // 3
    hist_ratio.at(0)->GetXaxis()->SetTitleSize(17); // 17
    hist_ratio.at(0)->GetXaxis()->SetTitleFont(46);
    hist_ratio.at(0)->GetYaxis()->SetNdivisions(5, 0, 0, kFALSE);
    hist_ratio.at(0)->GetYaxis()->SetRangeUser(0.3, 1.7);
    hist_ratio.at(0)->GetYaxis()->SetMaxDigits(1);
    hist_ratio.at(0)->GetYaxis()->SetTitle("Ratio to MC");
    hist_ratio.at(0)->GetYaxis()->SetTitleSize(13); // 13
    hist_ratio.at(0)->GetYaxis()->SetTitleFont(44);
    hist_ratio.at(0)->GetYaxis()->SetLabelSize(15); // new
    hist_ratio.at(0)->GetYaxis()->SetTitleOffset(2);
    hist_ratio.at(0)->SetTitle(" ");
    hist_ratio.at(0)->GetXaxis()->SetTitle(x_axis_name);
    hist_ratio.at(0)->GetXaxis()->SetRangeUser(0,7);

    // Draw the ratios on the ratio pad
    for (unsigned int y=0; y < hist.size(); y++ ) {
   
        hist_ratio.at(y)->Draw("hist,same");
    }

    h_error_hist_ratio->Draw("E2,same");
    h_error_hist_noDetvar_ratio->SetFillColorAlpha(kRed+2, 0.15);
    h_error_hist_noDetvar_ratio->Draw("E2,same");

    h_error_hist_data->Draw("PE, same");

    if (plotdata)
        _util.Draw_Data_POT(c, _util.config_v.at(_util.k_Run1_Data_POT), 0.45, 0.915, 0.45, 0.915);

    //---------------------------------------------------------------
    // draw final canvas as pdf

    c->Print("plots/WireModPaperPlot.pdf"); 

    // close the canvas to avoid warning messages on the terminal
    c->Close();  


    if (plotdata){
        delete h_data;
        delete h_dirt;
        delete h_ext;
    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::AddSysUncertainty(TH1D* h_error_hist, TH1D* h_ext, TH1D* h_dirt, std::string histname, std::string cut_name, std::string label, std::string mode){
        
    TH1D  *h_sys;

    TFile *file_sys_uncertainties = TFile::Open("files/run1_sys_var.root", "READ");

    // The error is on the MC events -- so comes from reweighting or detvar
    if (mode == "MC"){

        _util.GetHist(file_sys_uncertainties, h_sys, Form("%s/%s/%s", cut_name.c_str(), label.c_str(), histname.c_str()) );
        
        // loop over the bins in h_error_hist
        for (int i = 1; i <= h_error_hist->GetNbinsX() ; i++){

            double bin_error = h_error_hist->GetBinError(i);

            // Need to subtract the beam off and dirt
            double bin_content = h_error_hist->GetBinContent(i) - h_ext->GetBinContent(i) - h_dirt->GetBinContent(i);

            double sys_error = h_sys->GetBinContent(i) * bin_content;

            double tot_error = std::sqrt( bin_error*bin_error + sys_error*sys_error );

            h_error_hist->SetBinError(i, tot_error);

        }

        delete h_sys;
    }
    // Just add the uncertainty on the dirt -- e.g. add 100% uncertainty on the dirt
    else if (mode == "Dirt"){
        
        // loop over the bins in h_error_hist
        for (int i = 1; i <= h_error_hist->GetNbinsX() ; i++){

            double bin_error = h_error_hist->GetBinError(i);

            // Need to subtract the beam off and dirt
            double bin_content = h_dirt->GetBinContent(i);

            double sys_error = bin_content; // 100 % err on the dirt

            double tot_error = std::sqrt( bin_error*bin_error + sys_error*sys_error );

            h_error_hist->SetBinError(i, tot_error);

        }
    }
    // Add uncertainty on the total stack i.e. MC + dirt + EXT i.e.  in the case of POT counting
    else if (mode == "Stack"){

        // loop over the bins in h_error_hist
        for (int i = 1; i <= h_error_hist->GetNbinsX() ; i++){

            // std::cout <<histname  <<"  " << cut_name<< std::endl;

            double bin_error = h_error_hist->GetBinError(i);

            // Get the bin content
            double bin_content = h_error_hist->GetBinContent(i);

            double sys_error = 0.02 * bin_content; // 2 % err on the stack

            double tot_error = std::sqrt( bin_error*bin_error + sys_error*sys_error );

            h_error_hist->SetBinError(i, tot_error);

        }

    }
    else {
        std::cout << "Unknown systematics mode confugured!!" << std::endl;
    }

    file_sys_uncertainties->Close();

}
