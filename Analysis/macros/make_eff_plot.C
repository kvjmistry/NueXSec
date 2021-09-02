

// Declar functions
void IncreaseLabelSize(TH1D *h, TCanvas *c);
bool GetHist(TFile* f, TH1D* &h, TString string);
void MakeEfficiencyPlotByCutTot(std::string var_tot, std::string var_nue, std::string var_nuebar, std::string leg_tot, std::string leg_nue, std::string leg_nuebar, bool mask_title, bool mask_ax_label, const char* pri_ax_name, const char* printname);
void Draw_ubooneSim(TCanvas *c, double x1, double y1, double x2, double y2);


// MAIN
// -----------------------------------------------------------------------------
void make_eff_plot(){

    // Energy
    MakeEfficiencyPlotByCutTot("h_true_elec_E_rebin",     "h_true_elec_E_rebin_nue",     "h_true_elec_E_rebin_nuebar",     "#nu_{e} or #bar{#nu}_{e}", "#nu_{e}",          "#bar{#nu}_{e}",   true, false, "True E_{e} [GeV]; Efficiency",         "elec_E_rebin");
    
    // Angle
    MakeEfficiencyPlotByCutTot("h_eff_cosine_beta_rebin", "h_eff_cosine_beta_rebin_nue", "h_eff_cosine_beta_rebin_nuebar", "#nu_{e} or #bar{#nu}_{e}", "#nu_{e}",          "#bar{#nu}_{e}",   true, false, "True cos#beta_{e}; Efficiency",         "cosine_beta_rebin");

}
// -----------------------------------------------------------------------------


// Function to plot
void MakeEfficiencyPlotByCutTot(std::string var_tot, std::string var_nue, std::string var_nuebar, std::string leg_tot, std::string leg_nue, std::string leg_nuebar, bool mask_title, bool mask_ax_label, const char* pri_ax_name, const char* printname) {


    TFile *f_nuexsec = TFile::Open("nuexsec_run1_merged.root", "READ");

    // variables
    std::vector<std::string> cut_dirs = {"Unselected","dEdx_max_no_tracks"};
    enum enum_cut_dirs {k_unselected,k_dEdx_max_no_tracks,k_cuts_MAX}; 
    std::vector<std::string> cut_dirs_pretty = {"Unselected","dE/dx, 0 Tracks"};

    std::vector<TH1D *> hist_tot(k_cuts_MAX);    // The vector of histograms for nue plus nuebar / e- + e+
    std::vector<TH1D *> hist_nue(k_cuts_MAX);    // The vector of histograms for nue / e-
    std::vector<TH1D *> hist_nuebar(k_cuts_MAX); // The vector of histograms for nuebar / e+
    TH1D *h_clone_tot, *h_clone_nue, *h_clone_nuebar;

    // Helps determine what axes labels to draw 
    std::string var_string;

    // Loop over the classifications and get the histograms
    for (unsigned int i = 0; i < k_cuts_MAX; i++) {


        // MC
        GetHist(f_nuexsec, hist_tot.at(i),    Form("TEff/%s_%s", var_tot.c_str() ,   cut_dirs.at(i).c_str()));
        GetHist(f_nuexsec, hist_nue.at(i),    Form("TEff/%s_%s", var_nue.c_str() ,   cut_dirs.at(i).c_str()));
        GetHist(f_nuexsec, hist_nuebar.at(i), Form("TEff/%s_%s", var_nuebar.c_str() ,cut_dirs.at(i).c_str()));
        
        if (hist_tot.at(i) == NULL || hist_nue.at(i) == NULL || hist_nuebar.at(i) == NULL)
            return;
    }

    // Loop over the cuts and draw the efficiencies
    for (int p = 0; p < k_cuts_MAX; p++) {

        TCanvas * c = new TCanvas("c", "c", 500, 500);
        c->SetTopMargin(0.11);

        // Clone the histograms
        h_clone_tot    = (TH1D *) hist_tot.at(p)->Clone();
        h_clone_nue    = (TH1D *) hist_nue.at(p)->Clone();
        h_clone_nuebar = (TH1D *) hist_nuebar.at(p)->Clone();
        
        std::vector<double> eff_err_tot(h_clone_tot->GetNbinsX());
        std::vector<double> eff_err_nue(h_clone_nue->GetNbinsX());
        std::vector<double> eff_err_nuebar(h_clone_nuebar->GetNbinsX());
        
        // Get the bin errors based on binomial dist = sqrt(e/N*(1-e))) where e = n/N is the efficiency
        for (int bin = 0; bin < h_clone_tot->GetNbinsX(); bin++){
            
            // Get the errors for the nue plus nuebar hist
            double n = h_clone_tot->GetBinContent(bin+1)/0.09824; // selected = n
            double N = hist_tot.at(k_unselected)->GetBinContent(bin+1)/0.09824; // generated = N
            eff_err_tot.at(bin) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

            // nue
            n = h_clone_nue->GetBinContent(bin+1)/0.09824; // selected = n
            N = hist_nue.at(k_unselected)->GetBinContent(bin+1)/0.09824; // generated = N
            eff_err_nue.at(bin) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

            //nuebar
            n = h_clone_nuebar->GetBinContent(bin+1)/0.09824; // selected = n
            N = hist_nuebar.at(k_unselected)->GetBinContent(bin+1)/0.09824; // generated = N
            eff_err_nuebar.at(bin) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));

        }

        h_clone_tot->Divide(hist_tot.at(k_unselected));
        h_clone_nue->Divide(hist_nue.at(k_unselected));
        h_clone_nuebar->Divide(hist_nuebar.at(k_unselected));
        
        // Now set the bin errors after the divide
        for (int bin = 0; bin < h_clone_tot->GetNbinsX(); bin++){
            h_clone_tot->SetBinError(bin+1, eff_err_tot.at(bin));
            h_clone_nue->SetBinError(bin+1, eff_err_nue.at(bin));
            h_clone_nuebar->SetBinError(bin+1, eff_err_nuebar.at(bin));
        }
        
        
        h_clone_tot->SetStats(kFALSE);
        h_clone_nue->SetLineColor(kBlue+2);
        h_clone_nuebar->SetLineColor(kRed+2);
        h_clone_nue->SetLineWidth(3);
        h_clone_nuebar->SetLineWidth(3);
        h_clone_tot->SetTitle(Form("%s;%s", cut_dirs_pretty.at(p).c_str(), pri_ax_name));
        
        // Get rid of the ticks and the axes labels for single bin
        if (mask_ax_label) {
            h_clone_tot->GetXaxis()->SetLabelOffset(100);
            h_clone_tot->GetXaxis()->SetTickSize(0);
        }

        h_clone_tot->GetXaxis()->CenterTitle();
        h_clone_tot->GetYaxis()->SetRangeUser(0, 0.4);
        h_clone_tot->SetLineColor(kBlack);
        h_clone_tot->SetLineWidth(3);
        IncreaseLabelSize(h_clone_tot, c);
        c->SetLeftMargin(0.17);
        if (mask_title) h_clone_tot->SetTitle("");
        h_clone_tot->Draw("E same");
        h_clone_nue->Draw("E same");
        h_clone_nuebar->Draw("E same");

        std::size_t found = std::string(printname).find("multi"); // Look for "multi" in the name
        
        // Has multi in the name,so center the axes label
        if (found!=std::string::npos)
            h_clone_tot->GetXaxis()->CenterLabels();
        
        TLegend *leg = new TLegend(0.17, 0.75, 0.5, 0.89);
        leg->SetNColumns(1);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_clone_tot, leg_tot.c_str(), "le");  
        leg->AddEntry(h_clone_nue, leg_nue.c_str(), "le");  
        leg->AddEntry(h_clone_nuebar, leg_nuebar.c_str(), "le");  

        leg->Draw();

        TH1D* h_nue_clone = (TH1D*) hist_nue.at(p)->Clone();
        TH1D* h_nuebar_clone = (TH1D*) hist_nuebar.at(p)->Clone();
        
        found = std::string(printname).find("rebin"); // Look for "multi" in the name
        if (found!=std::string::npos){
            h_nue_clone->Scale(1.0, "width");
            h_nuebar_clone->Scale(1.0, "width");
        }

        
        h_nue_clone->SetLineWidth(2);
        h_nue_clone->SetLineStyle(2);

        
        h_nuebar_clone->SetLineWidth(2);
        h_nuebar_clone->SetLineStyle(3);
        h_nuebar_clone->SetLineColor(kRed+2);


        TLegend *leg2 = new TLegend(0.57, 0.75, 0.80, 0.85);
        leg2->SetNColumns(1);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->AddEntry(h_nue_clone, leg_nue.c_str(), "l");  
        leg2->AddEntry(h_nuebar_clone, leg_nuebar.c_str(), "l");  

        leg2->Draw();
        
        found = std::string(printname).find("elec_E"); // Look for "elec_E" in the name

        if (found!=std::string::npos){
            h_nue_clone->Scale(0.30 / (h_nue_clone->GetMaximum()));
            h_nuebar_clone->Scale(h_nue_clone->Integral() / (h_nuebar_clone->Integral()));
        }
        
        found = std::string(printname).find("cosine"); // Look for "cosine" in the name

        if (found!=std::string::npos){
            
            h_nuebar_clone->Scale(0.28 / (h_nuebar_clone->GetMaximum()));
            h_nue_clone->Scale(h_nuebar_clone->Integral() / (h_nue_clone->Integral()));
        }
        
        h_nue_clone->Draw("hist,same");

        h_nuebar_clone->Draw("hist,same");

        // Draw the run period on the plot
        Draw_ubooneSim(c, 0.73, 0.885, 0.73, 0.865);

        if (p !=0)
            c->Print(Form("TEff_%s_%s_combined.pdf", cut_dirs.at(p).c_str(), printname) );
        
        delete c;
        delete h_clone_tot;
        delete h_clone_nue;
        delete h_clone_nuebar;
    }

}
// -----------------------------------------------------------------------------
bool GetHist(TFile* f, TH1D* &h, TString string){
    f->cd();
    h = (TH1D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
void IncreaseLabelSize(TH1D *h, TCanvas *c){

    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);
}
// -----------------------------------------------------------------------------
void Draw_ubooneSim(TCanvas *c, double x1, double y1, double x2, double y2){
    c->cd();

    // 0.37, 0.92, 0.37, 0.92,

    TPaveText *pt;

    pt = new TPaveText(x1, y1, x2, y2,"NDC");
    pt->AddText("MicroBooNE Simulation");
    // pt->AddText("MicroBooNE Simulation In Progress");
    pt->SetTextColor(kBlack);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();

}
// -----------------------------------------------------------------------------