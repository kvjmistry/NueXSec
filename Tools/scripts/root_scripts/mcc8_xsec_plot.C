// Try improving the MCC8 x-sec plot
// Writing it in the MCC9 env so we can port this if we need it in the future.

void mcc8_xsec_plot(){

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c", 150, 350);
    gPad->SetLeftMargin(0.3);


    double x1[]  = {1.0};
    double y1[]  = {6.67e-39};
    double ex1[] = {0.0};
    double ey1[] = {0.267e-38};
    double ey2[] = {0.144e-38};
    
    // Data X-Sec with Stat Only
    TH1D* h_data = new TH1D("h_data", ";;#nu_{e} + #bar{#nu}_{e} CC Cross Section [cm^{2} / nucleon]", 1, 0, 1);
    
    // X-Axis
    h_data->GetXaxis()->SetRangeUser(0.0,1.0); 
    h_data->GetXaxis()->SetLabelOffset(999);
    h_data->GetXaxis()->SetLabelSize(0);
    h_data->GetXaxis()->SetTickLength(0);
    
    // Y-Axis
    h_data->GetYaxis()->SetRangeUser(0.2e-38,1.2e-38);
    h_data->GetYaxis()->CenterTitle();
    h_data->GetYaxis()->SetLabelSize(0.1);
    h_data->GetYaxis()->SetTitleSize(0.1);
   
    // h_data->SetLineWidth(2);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.7);
    h_data->SetLineColor(kBlack);

    TH1D * h_data_open = (TH1D*) h_data->Clone(); // for plotting the x-sec calculatiion without the open angle cut
    h_data_open->Fill(0.75, 7.40e-39);
    h_data_open->SetBinError(1, 1.52e-39);
    h_data_open->SetLineColor(kRed);
    h_data_open->SetMarkerColor(kRed);
    
    
    // Fill it
    h_data->Fill(0.5, 6.67e-39);
    h_data->SetBinError(1, 0.267e-38);
    h_data->Draw("E1,X0");

    // Systematic band
    TH1D * h_data_stat = (TH1D*) h_data->Clone();
    h_data_stat->SetBinError(1, 0.144e-38);
    h_data_stat->SetLineColor(kBlack);
    h_data_stat->Draw("E1,X0,same");


    

    // Genie nue + nuebar
    TH1D* h_genie_nue_nuebar = new TH1D("h_genie", "", 1, 0.0, 1.0);
    h_genie_nue_nuebar->Fill(0.5,6.9e-39 );
    h_genie_nue_nuebar->SetLineColor(kViolet-5);
    h_genie_nue_nuebar->SetLineWidth(3);
    h_genie_nue_nuebar->SetLineStyle(7);
    h_genie_nue_nuebar->Draw("hist,same");

    // Draw the open angle x-sec
    // h_data_open->Draw("E1,X0,same");

    // Draw the Legend
    TLegend *leg = new TLegend(0.35, 0.75, 0.70, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data, "#lower[-0.2]{Data (Stat. + Sys.)}",      "ep");
    // leg->AddEntry(h_data_open, "Data no O. A. (Stat.)",      "ep");
    // leg->AddEntry(g_xsec_sys,      "Stat. + Sys. Uncertainty",                  "f");
    // leg->AddEntry(h_data_stat,     "Data (Stat.)",                        "ep");
    // leg->AddEntry(g_xsec_nue,      "GENIE #nu_{e} CC Cross-Section",          "l");
    // leg->AddEntry(g_xsec_nuebar,   "GENIE #bar{#nu}_{e} CC Cross-Section",    "l");
    leg->AddEntry(h_genie_nue_nuebar,   "#lower[-0.2]{GENIE v2.12.2}",    "l");
    
    leg->Draw();

    gStyle->SetLegendTextSize(0.06);
   
    std::cout << gStyle->GetLegendTextSize() << std::endl;

    // pt->Draw();   


    // c->Range(0.0, 1.0, 0.2e-38,1.4e-38);
    // c->Update();

    TLatex *t = new TLatex(.40, .175, "#splitline{MicroBooNE NuMI}{Data 2.4#times10^{20} POT}");
    t->SetTextColor(kBlack);
    t->SetNDC();
    t->SetTextSize(1.5/30.);
    t->SetTextAlign(11);
    t->Draw();


    c->Print("plots/mcc8_nuexsec_generator_plot.pdf");


}