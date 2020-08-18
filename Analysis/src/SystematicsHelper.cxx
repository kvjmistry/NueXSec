#include "../include/SystematicsHelper.h"

// -----------------------------------------------------------------------------
void SystematicsHelper::Initialise(const char *_run_period, Utility _utility, const char* _mode){

    std::cout << "Initalising Systematics Helper..." << std::endl;
    _util = _utility;

    // Set the run period
    run_period = std::string(_run_period);

    // Off beam mode to compare bnb and numi off beam samples
    if (std::string(_mode) == "ext"){
        var_string = { "NuMI", "BNB" };
        mode = std::string(_mode);
    }

    // If we choose this mode then we actually want to use a different initialiser
    if (std::string(_mode) == "reweight"){
        InitialiseReweightingMode();
        return;
    }

    // Get the POT of the variations from the file
    GetPOT(_run_period);

    
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
    hist->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
    hist->GetYaxis()->SetTitleSize(12);
    hist->GetYaxis()->SetTitleFont(44);
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleOffset(1.5);
    hist->SetTitle(" ");

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
        if (mode == "default"){
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

            // Space Charge Corrected Y position comparision plot
            PlotVariationsEXT("h_reco_vtx_y_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_y_sce.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Y [cm]");


            // Space Charge Corrected X position comparision plot
            PlotVariationsEXT("h_reco_vtx_z_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_z_sce.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Z [cm]");
        
            // Leading Shower Phi
            PlotVariationsEXT("h_reco_leading_shower_phi", Form("plots/run%s/ext/comparisons/cuts/%s/reco_leading_shower_phi.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Phi [degrees]");

            // Leading Shower Theta
            PlotVariationsEXT("h_reco_leading_shower_theta", Form("plots/run%s/ext/comparisons/cuts/%s/reco_leading_shower_theta.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Theta [degrees]");


            // Shower Multiplicty
            PlotVariationsEXT("h_reco_shower_multiplicity", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shower_multiplicity.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Multiplicty");

            // Track Multiplicty
            PlotVariationsEXT("h_reco_track_multiplicity", Form("plots/run%s/ext/comparisons/cuts/%s/reco_track_multiplicity.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Track Multiplicity");

            
            // Topological Score
            PlotVariationsEXT("h_reco_topological_score", Form("plots/run%s/ext/comparisons/cuts/%s/reco_topological_score.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Topological Score");

            
            // Shower Track Fitter dedx Y plane
            PlotVariationsEXT("h_reco_shr_tkfit_dedx_y", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shr_tkfit_dedx_y.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Collection Plane dEdx (track fitter) [MeV/cm]");

            // Reco Electron Neutrino E
            PlotVariationsEXT("h_reco_nu_e", Form("plots/run%s/ext/comparisons/cuts/%s/reco_nu_e.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Neutrino Energy [GeV]");

            // Leading Shower Energy
            PlotVariationsEXT("h_reco_shower_energy_tot_cali", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shower_energy_tot_cali.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Leading Shower Energy [GeV]");

            // Flash PE
            PlotVariationsEXT("h_reco_flash_pe", Form("plots/run%s/ext/comparisons/cuts/%s/reco_flash_pe.pdf", run_period.c_str(), _util.cut_dirs.at(i).c_str()),
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
    // _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915, run_period);

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
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915, run_period);

    // Draw area normalisation
    Draw_Area_Norm(c);
    
    c->Print(print_name);



}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseReweightingMode(){

    gStyle->SetOptStat(0);

    // Load in the input file
    // Should we add more protection to this command??
    f_nuexsec = new TFile( Form("files/crosssec_run%s.root", run_period.c_str() ), "READ");

    InitialsePlotCV();

    for (unsigned int var = 0; var <  vars.size(); var++){

        CompareCVXSec(var);

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

        // Plot the unisims
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

        // Plot the multisims
        PlotReweightingModeMultisim("weightsGenie", var,  "GENIE", 500);
        PlotReweightingModeMultisim("weightsReint", var,  "Geant Reinteractions", 1000);
        PlotReweightingModeMultisim("weightsPPFX",  var,  "PPFX", 600);
        
    }
    

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeUnisim(std::string label, int var, std::string label_pretty){

    // Create the directory
    CreateDirectory("/Systematics/" + label + "/" + vars.at(var), run_period);

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
        _util.GetHist(f_nuexsec, h_universe.at(k_up).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), run_period.c_str(), label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, h_universe.at(k_dn).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), run_period.c_str(), label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
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
        
        h_universe.at(k_up).at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));
        h_universe.at(k_up).at(k)->GetXaxis()->SetTitle("");
        h_universe.at(k_up).at(k)->GetXaxis()->SetLabelSize(0);

        h_universe.at(k_up).at(k)->Draw("hist");
        cv_hist_vec.at(var).at(k)->Draw("hist,same");
        h_universe.at(k_dn).at(k)->Draw("hist,same");

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
        leg->AddEntry(h_universe.at(k_up).at(k), Form("%s +1 #sigma", label_pretty.c_str()), "l");
        leg->AddEntry(cv_hist_vec.at(var).at(k),         "CV", "l");
        leg->AddEntry(h_universe.at(k_dn).at(k), Form("%s -1 #sigma", label_pretty.c_str()), "l");
        leg->Draw();

        bottomPad->cd();
        
        // Up ratio to CV
        TH1D* h_err_up = (TH1D *)cv_hist_vec.at(var).at(k)->Clone("h_ratio_up");
        h_err_up->Add(h_universe.at(k_up).at(k), -1);
        h_err_up->Divide(cv_hist_vec.at(var).at(k));
        h_err_up->SetLineWidth(2);
        h_err_up->SetLineColor(kGreen+2);
        h_err_up->Scale(100);
        
        // Down ratio to CV
        TH1D* h_err_dn = (TH1D *)cv_hist_vec.at(var).at(k)->Clone("h_ratio_dn");
        h_err_dn->Add(h_universe.at(k_dn).at(k), -1);
        h_err_dn->Divide(cv_hist_vec.at(var).at(k));
        h_err_dn->SetLineWidth(2);
        h_err_dn->SetLineColor(kRed+2);
        h_err_dn->Scale(100);

        TH1D* h_err = (TH1D *)cv_hist_vec.at(var).at(k)->Clone("h_ratio");
        h_err->Divide(cv_hist_vec.at(var).at(k));

        SetRatioOptions(h_err_up);
        h_err_up->GetYaxis()->SetTitle("\% change from CV");
        h_err_up->Draw("hist,same, text00");
        h_err_up->SetMarkerSize(1.6);
        h_err_dn->Draw("hist,same, text00");
        h_err_dn->SetMarkerSize(1.6);
        
        
        gStyle->SetPaintTextFormat("4.2f");
        h_err->Draw("hist,same");

        c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", run_period.c_str(), label.c_str(), vars.at(var).c_str(), run_period.c_str(), label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_err;
        delete h_err_up;
        delete h_err_dn;
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeMultisim(std::string label, int var, std::string label_pretty, int universes){

    // Create the directory
    CreateDirectory("/Systematics/" + label + "/" + vars.at(var), run_period);

    std::vector<std::vector<TH1D*>> h_universe;
    std::vector<std::vector<TH1D*>> h_err;


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
            _util.GetHist(f_nuexsec, h_universe.at(uni).at(k), Form( "%s/%s/h_run%s_%s_%i_%s_%s", label.c_str(), vars.at(var).c_str(), run_period.c_str(), label.c_str(), uni ,vars.at(var).c_str(), xsec_types.at(k).c_str()));

            // Customise
            h_universe.at(uni).at(k)->SetLineWidth(1);
            h_universe.at(uni).at(k)->SetLineColor(kAzure+5);
            
        }
    }

    // Create vector of systematic error for each histogram
    std::vector<std::vector<double>> sys_err;
    sys_err.resize(xsec_types.size());
    
    for (unsigned int i = 0; i < sys_err.size(); i++){
        sys_err.at(i).resize(cv_hist_vec.at(var).at(0)->GetNbinsX(), 0.0);
    }

    // We now want to get the standard deviation of all universes wrt to the cv


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

        // Loop over universes
        for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
            h_universe.at(uni).at(k)->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            h_universe.at(uni).at(k)->Draw("hist,same");
            if (scale_val < h_universe.at(uni).at(k)->GetMaximum()) scale_val = h_universe.at(uni).at(k)->GetMaximum();
            h_universe.at(uni).at(k)->GetXaxis()->SetTitle("");
            h_universe.at(uni).at(k)->GetXaxis()->SetLabelSize(0);

        }

        if (xsec_types.at(k) == "mc_xsec") {
            // cv_hist_vec_clone.at(k)->SetFillColorAlpha(kRed+2, 0.00);
            // cv_hist_vec_clone.at(k)->SetFillColor(kBlack);
            // cv_hist_vec_clone.at(k)->SetFillStyle(3444);
            cv_hist_vec_clone.at(k)->SetLineColor(kBlack);
            cv_hist_vec_clone.at(k)->SetFillStyle(0);
            // cv_hist_vec_clone.at(k)->SetMarkerStyle(47);
            // cv_hist_vec_clone.at(k)->SetMarkerSize(0.4);
            cv_hist_vec_clone.at(k)->SetMarkerColor(kBlack);
            cv_hist_vec_clone.at(k)->Draw("E2,hist,same");
        }
        else { 
            cv_hist_vec_clone.at(k)->Draw("E,same");
        }

        

        
        

        h_universe.at(0).at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec){
            // h_universe.at(0).at(k)->GetYaxis()->SetRangeUser(0, 0.5e-39);
        }

        TLegend *leg = new TLegend(0.6, 0.65, 0.95, 0.8);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_universe.at(0).at(k), Form("%s - All Universes", label_pretty.c_str()), "l");
        if (xsec_types.at(k) == "mc_xsec")leg->AddEntry(cv_hist_vec_clone.at(k),           "Central Value", "ef");
        else leg->AddEntry(cv_hist_vec_clone.at(k),           "Central Value", "le");
        leg->Draw();

        bottomPad->cd();

        // Up percent diff to CV
        TH1D* h_err = (TH1D *)cv_hist_vec_clone.at(k)->Clone("h_err");
        
        // Loop over the bins in the up error, and set the bin content to be the percent difference
        for (int g = 1; g < h_err->GetNbinsX()+1; g++){
            h_err->SetBinContent(g, 100 * h_err->GetBinError(g)/h_err->GetBinContent(g));
        }
        h_err->SetLineWidth(2);
        h_err->SetLineColor(kAzure+5);
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
        
        h_err->GetYaxis()->SetTitle("\% Uncertainty");
        h_err->Draw("hist, text00");
        h_err->SetMarkerSize(3);
        gStyle->SetPaintTextFormat("4.2f");

        c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", run_period.c_str(), label.c_str(), vars.at(var).c_str(),  run_period.c_str(), label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;

    }

        
    
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareCVXSec(int var){

    
    TH1D* h_dataxsec = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_temp");
    TH1D* h_mcxsec   = (TH1D*) cv_hist_vec.at(var).at(k_xsec_mcxsec)  ->Clone("h_data_xsec_temp");


    TPad *topPad;
    TPad *bottomPad;

    TCanvas * c = new TCanvas("c", "c", 500, 500);
    topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    _util.SetTPadOptions(topPad, bottomPad);

    h_dataxsec->SetLineColor(kGreen+2);
    h_mcxsec  ->SetLineColor(kRed+2);

    // h_dataxsec->GetYaxis()->SetRangeUser(0, 0.5e-39);
    if (vars.at(var) == "integrated") h_dataxsec->GetYaxis()->SetRangeUser(0, 10e-39);

    h_dataxsec->GetYaxis()->SetLabelSize(0.04);
    h_dataxsec->GetYaxis()->SetTitleSize(14);
    h_dataxsec->GetYaxis()->SetTitleFont(44);
    h_dataxsec->GetYaxis()->SetTitleOffset(1.5);
    h_dataxsec->GetXaxis()->SetTitle("");
    h_dataxsec->GetXaxis()->SetLabelSize(0);

    h_dataxsec->Draw("hist,E");
    h_mcxsec->Draw("hist,same");
    TH1D* h_mcxsec_clone = (TH1D *)h_mcxsec->Clone("h_mc_clone");
    h_mcxsec_clone->SetFillColorAlpha(12, 0.15);
    h_mcxsec_clone->Draw("E2,same");
    

    TLegend *leg = new TLegend(0.6, 0.7, 0.95, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_dataxsec, "Data CV", "le");
    leg->AddEntry(h_mcxsec_clone,   "MC CV", "lf");
    leg->Draw();


    bottomPad->cd();
        
    // The percent difference of mc wrt data
    TH1D* h_err = (TH1D *)h_dataxsec->Clone("h_ratio");
    h_err->Add(h_mcxsec, -1);
    h_err->Divide(h_dataxsec);
    h_err->SetLineWidth(2);
    h_err->SetLineColor(kGreen+2);
    h_err->Scale(100);

    bottomPad->SetGridy(kFALSE);
    

    SetRatioOptions(h_err);
    h_err->SetLineColor(kBlack);
    h_err->GetYaxis()->SetTitleSize(11);
    h_err->GetYaxis()->SetRangeUser(-100, 100);
    h_err->GetYaxis()->SetTitle("Data - MC / Data [\%]");
    h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
    h_err->Draw("hist,same");

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915, run_period);


    c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison.pdf", run_period.c_str(), vars.at(var).c_str(), run_period.c_str(), vars.at(var).c_str() ));
    delete c;
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
            _util.GetHist(f_nuexsec, cv_hist_vec.at(var).at(k), Form( "CV/%s/h_run%s_CV_0_%s_%s", vars.at(var).c_str(), run_period.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

            if (cv_hist_vec.at(var).at(k) == NULL) std::cout << "Failed to get the histogram!" << std::endl;

            // Customise
            cv_hist_vec.at(var).at(k)->SetLineWidth(2);
            cv_hist_vec.at(var).at(k)->SetLineColor(kBlack);
        }

        // Create the CV directory and draw the CV
        CreateDirectory("/Systematics/CV/" + vars.at(var) + "/", run_period);
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
            _util.Draw_Run_Period(c_cv, 0.86, 0.92, 0.86, 0.92, run_period);

            c_cv->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_%s.pdf", run_period.c_str(), vars.at(var).c_str(), run_period.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

            delete c_cv;
            delete leg;
        }
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareVariationXSec(std::string label, int var, std::string label_pretty){

    
        // Create the directory
    CreateDirectory("/Systematics/" + label + "/" + vars.at(var), run_period);

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
        _util.GetHist(f_nuexsec, h_universe.at(k_up).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), run_period.c_str(), label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        if (k == k_xsec_mcxsec) h_universe.at(k_up).at(k)->SetLineStyle(7);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, h_universe.at(k_dn).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), run_period.c_str(), label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

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
    h_err_up->GetYaxis()->SetRangeUser(-100, 100);
    h_err_up->GetYaxis()->SetTitle("\% change of Data to MC");
    h_err_up->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
    h_err_up->Draw("hist,same");
    h_err_dn->Draw("hist,same");
    h_err->Draw("hist,same");

    c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_data_mc_comparison.pdf", run_period.c_str(), label.c_str(), vars.at(var).c_str(), run_period.c_str(), label.c_str(), vars.at(var).c_str() ));

    delete c;
}