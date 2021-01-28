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
    h_data->GetYaxis()->SetRangeUser(0.22e-38,1.22e-38);
    h_data->GetYaxis()->CenterTitle();
    h_data->GetYaxis()->SetLabelSize(0.1);
    h_data->GetYaxis()->SetTitleSize(0.1);
   
    // h_data->SetLineWidth(2);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.7);
    h_data->SetLineColor(kBlack);
    
    // Fill it
    // h_data->Fill(0.5, 6.67e-39);
    // h_data->SetBinError(1, 0.267e-38);
    h_data->Fill(0.5, 6.8426915e-39); // new in FV flux
    h_data->SetBinError(1, 6.8426915e-39 * 0.40); // new in FV flux
    h_data->Draw("E1,X0");

    // Statistical band
    TH1D * h_data_stat = (TH1D*) h_data->Clone();
    // h_data_stat->SetBinError(1, 0.144e-38);
    h_data_stat->SetBinError(1, 6.8426915e-39 * 0.22 ); // new in FV flux
    h_data_stat->SetLineColor(kBlack);
    h_data_stat->Draw("E1,X0,same");


    // Genie v12.2.2 nue + nuebar
    TH1D* h_genie_v2_nue_nuebar = new TH1D("h_genie_v2", "", 1, 0.0, 1.0);
    // h_genie_v2_nue_nuebar->Fill(0.5,7.19925e-39 );
    h_genie_v2_nue_nuebar->Fill(0.5,7.3125100e-39 ); // with FV flux
    h_genie_v2_nue_nuebar->SetLineColor(kViolet-5);
    h_genie_v2_nue_nuebar->SetLineWidth(3); 
    h_genie_v2_nue_nuebar->SetLineStyle(7);
    h_genie_v2_nue_nuebar->Draw("hist,same");

    // Genie v3 nue + nuebar
    TH1D* h_genie_v3_nue_nuebar = new TH1D("h_genie_v3", "", 1, 0.0, 1.0);
    // h_genie_v3_nue_nuebar->Fill(0.5,5.5228738e-39 );
    h_genie_v3_nue_nuebar->Fill(0.5,5.5711475e-39 ); // with FV flux
    h_genie_v3_nue_nuebar->SetLineColor(kBlue+2);
    h_genie_v3_nue_nuebar->SetLineWidth(3);
    h_genie_v3_nue_nuebar->SetLineStyle(8);
    h_genie_v3_nue_nuebar->Draw("hist,same");

    // NuWro nue + nuebar
    TH1D* h_genie_NuWro_nue_nuebar = new TH1D("h_nuwro_v2", "", 1, 0.0, 1.0);
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
    leg->AddEntry(h_genie_v2_nue_nuebar,   "GENIE v2.12.2",    "l");
    leg->AddEntry(h_genie_v3_nue_nuebar,   "GENIE v3.0.6",    "l");
    leg->AddEntry(h_genie_NuWro_nue_nuebar,   "NuWro v19.02.1",    "l");
    
    leg->Draw();

    gStyle->SetLegendTextSize(0.06);
   
    std::cout << gStyle->GetLegendTextSize() << std::endl;

    TLatex *t = new TLatex(.34, .145, "#splitline{MicroBooNE NuMI}{Data 2.4#times10^{20} POT}");
    t->SetTextColor(kBlack);
    t->SetNDC();
    t->SetTextSize(2.0/30.);
    t->SetTextAlign(11);
    t->Draw();


    c->Print("plots/mcc8_nuexsec_generator_plot.pdf");


}