// Get Histograms
void GetHist(TFile* f, TH1D* &h, TString string);

// Define a tuple for the variations

// var name, var name pretty, colour
std::vector<std::tuple<std::string, std::string, int>> tuple_label = {
    std::make_tuple("CV", "CV", 1),
    std::make_tuple("WireModX", "WM X", 2),
    std::make_tuple("WireModYZ", "WM YZ", 4),
    std::make_tuple("WireModThetaXZ", "WM Theta XZ", 6),
    std::make_tuple("WireModThetaYZ_withoutSigmaSplines", "WM Theta YZ", 8)
};


// Recreate the WireMod Paper plot
void WireModPaperPlot(){

    gStyle->SetOptStat(0);

    // Load in the TFile
    TFile *f_WM = TFile::Open("WireModPaperPlot.root", "READ");

    // Top Pad Histograms
    
    std::vector<TH1D*> hist(tuple_label.size()); // WM var hists
    TH1D *h_data;                                // beam on data
    TH1D *h_error_hist;                          // Total Uncertainty (grey band)

    // Bottom Pad Histograms
    std::vector<TH1D*> hist_ratio(tuple_label.size()); // WM var ratio hists
    TH1D *h_data_ratio;                                // beam on data

    TH1D* h_error_hist_ratio;          // The total error band in ratio
    TH1D* h_error_hist_noDetvar_ratio; // The total error band without WM


    // -- Get the Histograms---

    // Get the histograms and Customise
    for (unsigned int y=0; y < hist.size(); y++ ) {
        GetHist(f_WM, hist.at(y), Form("%s",  std::get<0>(tuple_label.at(y)).c_str()));
        GetHist(f_WM, hist_ratio.at(y), Form("%s_ratio",  std::get<0>(tuple_label.at(y)).c_str()));

        hist.at(y)->SetLineWidth(2);
        hist_ratio.at(y)->SetLineWidth(2);
        hist.at(y)->SetLineColor(std::get<2>(tuple_label.at(y)));
        hist_ratio.at(y)->SetLineColor(std::get<2>(tuple_label.at(y)));
    }

    GetHist(f_WM, h_data, "h_data");
    GetHist(f_WM, h_error_hist, "h_error_hist");
    GetHist(f_WM, h_data_ratio, "h_error_hist_ratio_data");
    GetHist(f_WM, h_error_hist_ratio, "h_error_hist_ratio");
    
    
    GetHist(f_WM, h_error_hist_noDetvar_ratio, "h_error_hist_ratio_noDetvar");
    h_error_hist_noDetvar_ratio->SetFillColorAlpha(kRed+2, 0.15);

    //  --- Plotting ---

    // Create the Canvas
    TCanvas * c      = new TCanvas("","",500,500);
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    topPad->SetBottomMargin(0.05);
    topPad->SetTopMargin(0.15);
    bottomPad->SetTopMargin(0.04);
    bottomPad->SetBottomMargin(0.3);
    bottomPad->SetGridy();
    topPad->SetLeftMargin(0.15);
    topPad->SetRightMargin(0.1);
    bottomPad->SetLeftMargin(0.15);
    bottomPad->SetRightMargin(0.1);
    topPad->Draw();
    bottomPad->Draw();
    topPad->cd();

    // Legend
    // on top of the topCanvs to avoid overlaping the plot
    TLegend *leg = new TLegend(0.1686747,0.7233083,0.8795181,0.8406015,NULL,"brNDC");
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // variable used to track the max bin content to later scale the plot to avoid cutting information
    
    // Set the histogram scales and add to legend
    for (unsigned int y=0; y < hist.size(); y++ ) {
        
        // CV
        if (y == 0) {
            leg->AddEntry(h_error_hist, "CV (Tot Unc.)", "lf");        
        }
        else {
            hist.at(y)->SetLineColor(std::get<2>(tuple_label.at(y)));  // change color of the histogram
            hist_ratio.at(y)->SetLineColor(std::get<2>(tuple_label.at(y)));
            leg->AddEntry(hist.at(y),std::get<1>(tuple_label.at(y)).c_str(), "l"); // add histogram to legend
            
        }
    }

    // Add the beam-on to the legend
    leg->AddEntry(h_data, "Beam-On", "lep"); // add histogram to legend

    // -----------------------------------------------------------------
    // Drawing histograms on top pad
    for (unsigned int y=0; y < hist.size(); y++ ) {
        hist.at(0)->GetXaxis()->SetRangeUser(0,7);
        hist.at(y)->Draw("hist, same");
    }
    
    // -----------------------------------------------------------------
    
    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("E2, same");

    // drawing CV again to make sure it is on top of everything else
    hist.at(0)->Draw("hist, same");
    h_data->Draw("same PE");
    
    // drawing legend
    leg->Draw();

    // -----------------------------------------------------------------
    // Drawing the total detector systematic uncertainty on the bottom pad

    bottomPad->cd();

    // Draw the ratios on the ratio pad
    for (unsigned int y=0; y < hist.size(); y++ ) {
   
        hist_ratio.at(y)->Draw("hist,same");
    }

    h_error_hist_ratio->SetFillColorAlpha(12, 0.15);
    h_error_hist_ratio->Draw("E2,same");
    h_error_hist_noDetvar_ratio->SetFillColorAlpha(kRed+2, 0.15);
    h_error_hist_noDetvar_ratio->Draw("E2,same");

    h_data_ratio->Draw("PE, same");

    
    // Draw the POT on the plot
    c->cd();
    double POT = 2.0e20;
    TPaveText * pt= new TPaveText( 0.45, 0.915, 0.45, 0.915, "NDC");
    pt->AddText("MicroBooNE NuMI Data: 2.0#times10^{20} POT");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();

    c->Print("WireModPaperPlot.pdf"); 
    c->Update();

}


void GetHist(TFile* f, TH1D* &h, TString string){
    f->cd();
    h = (TH1D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return;
    }
    else {
        return;
    }
}