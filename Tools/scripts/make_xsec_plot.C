// Script to make the cross section plot
#include "functions.h"

void make_xsec_plot(){

    gStyle->SetOptStat(0);

    TFile *f;

    TH1D *h_nue, *h_nuebar;
    TGraph *genieXsecNueCC;
    TGraph *genieXsecNuebarCC;

    GetFile(f, "/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root");
    
    // Get the Flux Histograms
    GetHist(f, h_nue,    "nueFluxHisto");
    GetHist(f, h_nuebar, "anueFluxHisto");

    const char* genieXsecPath = "/cvmfs/uboone.opensciencegrid.org/products/genie_xsec/v2_12_0/NULL/DefaultPlusMECWithNC/data";
    
    if ( !genieXsecPath ) {
        std::cout << "$(GENIEXSECPATH) not defined." << std::endl;
        std::cout << "Please setup *genie_xsec*." << std::endl; 
    }
    else {
        TString genieXsecFileName = genieXsecPath;
        genieXsecFileName += "/xsec_graphs.root";
        TFile *genieXsecFile = new TFile(genieXsecFileName,"READ");
        genieXsecNueCC     = (TGraph *) genieXsecFile->Get("nu_e_Ar40/tot_cc");
        genieXsecNuebarCC  = (TGraph *) genieXsecFile->Get("nu_e_bar_Ar40/tot_cc");
        
        genieXsecFile->Close();
    }

    int num_bins = genieXsecNueCC->GetN()+1;
    std::cout << "number of bins: " << num_bins<< std::endl;

    TH1D *h_spline_nue    = new TH1D("h_spline_nue", ";Energy [GeV]; Cross-Section",num_bins, 0, 125);
    TH1D *h_spline_nuebar = new TH1D("h_spline_nuebar", ";Energy [GeV]; Cross-Section",num_bins, 0, 125);

    // Convert TGraph to histogram for nue and nuebar
    
    // Nue
    for (int i=0; i< genieXsecNueCC->GetN()+0; i++){
        double x,y;
        genieXsecNueCC->GetPoint(i, x, y);
        h_spline_nue->Fill(x, y/40.0);
    }

    // Nuebar
    for (int i=0; i< genieXsecNuebarCC->GetN(); i++){
        double x,y;
        genieXsecNuebarCC->GetPoint(i, x, y);
        h_spline_nuebar->Fill(x, y/40.0);
    }


    TCanvas *c = new TCanvas();
    
    // Nue flux
    h_nue->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu}_{e} / cm^{2} / 6 #times 10^{20} POT");
    h_nue->GetXaxis()->CenterTitle();
    h_nue->GetYaxis()->CenterTitle();

    h_nue->GetXaxis()->SetLabelFont(42);
    h_nue->GetXaxis()->SetLabelSize(0.04);
    h_nue->GetXaxis()->SetTitleSize(0.04);
    h_nue->GetYaxis()->SetLabelFont(42);
    h_nue->GetYaxis()->SetLabelSize(0.04);
    h_nue->GetYaxis()->SetTitleSize(0.04);

    h_nue->GetXaxis()->SetRangeUser(0,3);
    h_nue->GetYaxis()->SetRangeUser(0,150.0e6);
    h_nue->SetLineColor(kBlue+2);
    h_nue->SetFillColor(17);



    h_nue->GetXaxis()->SetTitleSize(17);
    h_nue->GetXaxis()->SetTitleFont(46);

    h_nue->GetXaxis()->SetTitleSize(18);
    h_nue->GetXaxis()->SetTitleFont(46);

    h_nue->Draw("hist");
    
    // Nuebar flux
    h_nuebar->GetXaxis()->SetRangeUser(0,3);
    h_nuebar->SetLineColor(kGreen+2);
    h_nuebar->SetFillColor(16);
    h_nuebar->Draw("hist,same");
    

    // Now setup the twin axis
    gPad->SetRightMargin(0.17 );

    Float_t rightmax =  3.0;
    Float_t scale = 150.0e6/rightmax;
    
    // nue spline
    h_spline_nue->SetLineWidth(2);
    h_spline_nue->SetLineColor(kBlue+2);
    h_spline_nue->SetLineStyle(3);
    h_spline_nue->Scale(scale);
    h_spline_nue->Draw("hist,L,P,same");

    //nuebar spline
    h_spline_nuebar->SetLineWidth(2);
    h_spline_nuebar->SetLineColor(kGreen+2);
    h_spline_nuebar->SetLineStyle(3);
    h_spline_nuebar->Scale(scale);
    h_spline_nuebar->Draw("hist,L,P,same");

    // The second axis
    TGaxis *axis = new TGaxis(3.0, 0, 3.0, 150.0e6, 0, rightmax, 510, "+L");
    axis->SetTitle("#nu_{e} + #bar{#nu}_{e} CC Cross-Section [10^{-38} cm^{2}]");
    axis->SetTitleOffset(1.1);
    axis->SetLineColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->SetTitleColor(kBlack);
    axis->SetTextFont(42);
    axis->SetLabelFont(42);
    axis->CenterTitle();

    axis->Draw();

    // Get the energy values
    double average_num = 0; // flux numerator
    double average_den = h_nue->Integral() +  h_nuebar->Integral(); // flux denominator

    int nue_flux_bins  = h_nue->GetNbinsX();
    int anue_flux_bins = h_nuebar->GetNbinsX();

    for (int k = 0; k < nue_flux_bins; k++ ){
        average_num = average_num + (h_nue->GetBinContent(k) * h_nue->GetBinCenter(k));
    }

    for (int k = 0; k < anue_flux_bins; k++ ){
        average_num = average_num + (h_nuebar->GetBinContent(k) * h_nuebar->GetBinCenter(k));
    }

    double average = average_num / average_den;
    std::cout << "Average Energy: " << average << std::endl;

    // Sum the fluxes
    TH1D *summed_flux = new TH1D("summed_flux", "summed_flux", 400, 0, 20);
    summed_flux->Add(h_nue, 1);
    summed_flux->Add(h_nuebar, 1);

    
    double summed_flux_integral = summed_flux->Integral();
    
    std::cout <<  "summed_flux_integral: " << summed_flux_integral << std::endl;

    double average_bin = 0;
    for (int bin = 0 ; bin < nue_flux_bins; bin++ ){
        
        if (summed_flux->GetBinCenter(bin) <= average){
            continue;
        }

        // here we should get the bin with the average in it
        average_bin = bin;
        break;
    }

    double max_sum = 0;
    double max_bin_val = 0;
    for (int bin = average_bin; bin < nue_flux_bins; bin++){
        
        max_sum = max_sum + summed_flux->GetBinContent(bin);
        
        if (max_sum / summed_flux_integral >= 0.34){
            max_bin_val = summed_flux->GetBinCenter(bin);
            break;
        }
    }
       

    double min_sum = 0;
    double min_bin_val = 0;
    for (int bin = average_bin; bin > 0; bin--){
        
        // print "Bin: ", bin
        min_sum = min_sum + summed_flux->GetBinContent(bin);
        
        if (min_sum / summed_flux_integral >= 0.34){
            min_bin_val = summed_flux->GetBinCenter(bin);
            break;
        }
    }

    std::cout <<"max_bin_val: " << max_bin_val << std::endl;
    std::cout <<"min_bin_val: " << min_bin_val << std::endl;

    double energy_err_min = average - min_bin_val;
    double energy_err_max = max_bin_val - average;


    // Now draw the data cross section
    double data_xsec = 0.467 * scale; // Scale is to take it to the flux axes

    double data_stat_err_plus = 0.101 * scale;
    double data_stat_err_minus = 0.101 * scale;

    double data_stat_plus_sys_err_plus = 0.1907 * scale;
    double data_stat__plus_sys_err_minus = 0.1907 * scale;

    const Int_t n = 2;
    Double_t x[n]   = {average, average};
    Double_t y[n]   = {data_xsec,data_xsec};
    Double_t exl[n] = {energy_err_min, energy_err_min};
    Double_t eyl[n] = {data_stat_err_minus, data_stat__plus_sys_err_minus};
    Double_t exh[n] = {energy_err_max, energy_err_max};
    Double_t eyh[n] = {data_stat_err_plus, data_stat_plus_sys_err_plus};
    
    auto gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
    
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(kFullCircle);
    gr->SetLineWidth(2);
    gr->Draw("L,P,same");

    // Draw the Legend
    TLegend *leg = new TLegend(0.17, 0.5, 0.48, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_nue,             "NuMI #nu_{e} Flux",                       "l");
    leg->AddEntry(h_nuebar,          "NuMI #bar{#nu}_{e} Flux",                 "l");
    leg->AddEntry(h_spline_nue,      "GENIE #nu_{e} CC Cross-Section",          "l");
    leg->AddEntry(h_spline_nuebar,   "GENIE #bar{#nu}_{e} CC Cross-Section",    "l");
    leg->AddEntry(gr, "Data #sigma_{#nu_{e} + #bar{#nu}_{e}} (Stat. + Sys.)",   "lep");
    leg->Draw();

    TPaveText *pt;

    pt = new TPaveText(0.27, 0.92, 0.27, 0.92,"NDC");
    pt->AddText("MicroBooNE");
    pt->SetTextColor(kBlack);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.04);
    pt->Draw();


    c->Print("../../Analysis/plots/flux_combo_morebins_v3.pdf");



    // Now create the other flux integrated cross sec plot
    TCanvas *c2  = new TCanvas();
    double x1[]  = {0, 1};
    double y1[]  = {4.65e-39, 4.65e-39};
    double ex1[] = {0.0, 0.0};
    double ey1[] = {0.1907e-38, 0.1907e-38};
    double ey2[] = {0.101e-38, 0.101e-38};
    
    // Stats band
    auto g_xsec_sys = new TGraphErrors(2, x1, y1, ex1, ey1);
    g_xsec_sys->SetFillColor(kGray);
    g_xsec_sys->SetFillColorAlpha(12, 0.15);
    g_xsec_sys->SetTitle(";;#nu_{e} + #bar{#nu}_{e} CC Cross-Section [cm^{2}]");
    g_xsec_sys->GetXaxis()->SetRangeUser(0,1);
    g_xsec_sys->GetYaxis()->SetRangeUser(0,1.0e-38);
    g_xsec_sys->GetXaxis()->SetLabelOffset(10);
    gStyle->SetTickLength(0.00,"x"); 
    g_xsec_sys->GetXaxis()->SetLabelOffset(999);
    g_xsec_sys->GetXaxis()->SetLabelSize(0);
    g_xsec_sys->GetXaxis()->SetTickLength(0);
    g_xsec_sys->GetYaxis()->CenterTitle();
    g_xsec_sys->GetYaxis()->SetLabelSize(0.04);
    g_xsec_sys->GetYaxis()->SetTitleSize(0.04);
    g_xsec_sys->Draw("a3");


    // Systematic band
    auto g_xsec_stat = new TGraphErrors(2, x1, y1, ex1, ey2);
    g_xsec_stat->SetFillColorAlpha(46, 0.15);
    g_xsec_stat->SetTitle(";;#nu_{e} + #bar{#nu}_{e} CC Cross-Section [10^{-38} cm^{2}]");
    g_xsec_stat->GetXaxis()->SetRangeUser(0,1);
    g_xsec_stat->GetYaxis()->SetRangeUser(0,1.0e-38);
    g_xsec_stat->Draw("3, same");

    // data xsec
    double y_xsec[]  = {4.65e-39, 4.65e-39};
    auto g_xsec = new TGraphErrors(2, x1, y_xsec);
    g_xsec->SetLineColor(kBlack);
    g_xsec->SetLineWidth(2);
    g_xsec->Draw("same");

    // Genie nue
    double y_nue[]  = {6.34569e-39, 6.34569e-39};
    auto g_xsec_nue = new TGraphErrors(2, x1, y_nue);
    g_xsec_nue->SetLineColor(kBlue+2);
    g_xsec_nue->SetLineWidth(2);
    g_xsec_nue->SetLineStyle(3);
    g_xsec_nue->Draw("same");

    // Genie nuebar
    double y_nuebar[]  = {2.24685e-39, 2.24685e-39};
    auto g_xsec_nuebar = new TGraphErrors(2, x1, y_nuebar);
    g_xsec_nuebar->SetLineColor(kGreen+2);
    g_xsec_nuebar->SetLineWidth(2);
    g_xsec_nuebar->SetLineStyle(3);
    g_xsec_nuebar->Draw("same");

    // Genie nue + nuebar
    double y_nue_nuebar[]  = {4.83e-39, 4.83e-39};
    auto g_xsec_nue_nuebar = new TGraphErrors(2, x1, y_nue_nuebar);
    g_xsec_nue_nuebar->SetLineColor(kViolet-5);
    g_xsec_nue_nuebar->SetLineWidth(2);
    g_xsec_nue_nuebar->SetLineStyle(7);
    g_xsec_nue_nuebar->Draw("same");

    // Draw the Legend
    TLegend *leg2 = new TLegend(0.57, 0.65, 0.88, 0.9);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(g_xsec, "Data #nu_{e} + #bar{#nu}_{e} CC Cross-Section",      "l");
    leg2->AddEntry(g_xsec_sys,      "Stat. + Sys. Uncertainty",                  "f");
    leg2->AddEntry(g_xsec_stat,     "Stat Uncertainty",                        "f");
    leg2->AddEntry(g_xsec_nue,      "GENIE #nu_{e} CC Cross-Section",          "l");
    leg2->AddEntry(g_xsec_nuebar,   "GENIE #bar{#nu}_{e} CC Cross-Section",    "l");
    leg2->AddEntry(g_xsec_nue_nuebar,   "GENIE #nu_{e} + #bar{#nu}_{e} CC Cross-Section",    "l");
    
    leg2->Draw();
   
    pt->Draw();   

    c2->Print("../../Analysis/plots/flux_combo_2_with_data_v3.pdf");

}
