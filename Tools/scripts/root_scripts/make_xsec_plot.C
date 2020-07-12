// Script to make the cross section plot
#include "functions.h"
#include "TGaxis.h"

void make_xsec_plot(){

    // This changes the plot to average splines/flux or not
    bool draw_averge = true;

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
    TH1D *h_spline_average = new TH1D("h_spline_average", ";Energy [GeV]; Cross-Section",num_bins, 0, 125);

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

    h_spline_average->Add(h_spline_nue,1);
    h_spline_average->Add(h_spline_nuebar,1);
    h_spline_average->Scale(0.5);


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

    if (!draw_averge) h_nue->Draw("hist");

    TH1D* h_nue_clone = (TH1D*)h_nue->Clone("h_nue_clone");
    
    // Nuebar flux
    h_nuebar->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu}_{e} / cm^{2} / 6 #times 10^{20} POT");
    h_nuebar->GetXaxis()->SetRangeUser(0,3);
    h_nuebar->SetLineColor(kGreen+2);
    h_nuebar->SetFillColor(16);
    if (!draw_averge)h_nuebar->Draw("hist,same");

    // Define the threshold 
    double threshold_energy =  0.250; // 0.75*0.2065
    std::cout << "Theshold Energy: " << threshold_energy*1000 << " MeV" << std::endl;
    
    double xbin_th = h_nue->GetXaxis()->FindBin(threshold_energy); // find the x bin to integrate from (threshold)
    double kdar_max = h_nue->GetXaxis()->FindBin( 0.225); // end of KDAR spectrum
    double flux_int_thresh_nue = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1);
    double flux_int_thresh_nuebar = h_nuebar->Integral( xbin_th, h_nuebar->GetNbinsX()+1);

    // Sum the fluxes
    TH1D *summed_flux = (TH1D*)h_nue->Clone("h_summed_flux");
    
    // summed_flux->Add(h_nue, 1);
    summed_flux->Add(h_nuebar, 1);

    TH1D* h_summed_flux_clone = (TH1D*)summed_flux->Clone("h_summed_flux_clone");
    // Here loop over the cloned histogram of nue, and zero out the flux after the threshold, this is so we can shade out the thresholded area
    for (unsigned int p=xbin_th; p < h_summed_flux_clone->GetNbinsX()+1; p++){
        h_summed_flux_clone->SetBinContent(p, 0);
    }

    summed_flux->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e} + #bar{#nu}_{e} / cm^{2} / 6 #times 10^{20} POT");
    summed_flux->GetXaxis()->CenterTitle();
    summed_flux->GetYaxis()->CenterTitle();
    summed_flux->GetXaxis()->SetLabelFont(42);
    summed_flux->GetXaxis()->SetLabelSize(0.04);
    // summed_flux->GetXaxis()->SetTitleSize(0.04);
    summed_flux->GetYaxis()->SetLabelFont(42);
    summed_flux->GetYaxis()->SetLabelSize(0.04);
    summed_flux->GetYaxis()->SetTitleSize(0.04);
    summed_flux->GetXaxis()->SetRangeUser(0,3);
    summed_flux->GetYaxis()->SetRangeUser(0,170.0e6);
    summed_flux->SetLineColor(kRed+2);
    summed_flux->SetFillColor(17);
    if (draw_averge) summed_flux->Draw("hist");
    
    

    // Now setup the twin axis
    gPad->SetRightMargin(0.17 );

    c->Update();

    Float_t rightmax =  3.0;
    Float_t scale =  gPad->GetUymax()/rightmax;
    
    // nue spline
    h_spline_nue->SetLineWidth(2);
    h_spline_nue->SetLineColor(kBlue+2);
    h_spline_nue->SetLineStyle(3);
    h_spline_nue->Scale(scale);
    if (!draw_averge)h_spline_nue->Draw("hist,L,P,same");

    //nuebar spline
    h_spline_nuebar->SetLineWidth(2);
    h_spline_nuebar->SetLineColor(kGreen+2);
    h_spline_nuebar->SetLineStyle(3);
    h_spline_nuebar->Scale(scale);
    if (!draw_averge)h_spline_nuebar->Draw("hist,L,P,same");

    // Averge spline
    h_spline_average->SetLineWidth(2);
    h_spline_average->SetLineColor(kRed+2);
    h_spline_average->SetLineStyle(3);
    h_spline_average->Scale(scale);
    h_spline_average->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e} + #bar{#nu}_{e} / cm^{2} / 6 #times 10^{20} POT");
    if (draw_averge) h_spline_average->Draw("hist,L,P,same");

    // The second axis
    TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax, 510, "+L");
    axis->SetTitle("#nu_{e} + #bar{#nu}_{e} CC Cross-Section [10^{-38} cm^{2}]");
    axis->SetTitleOffset(1.1);
    axis->SetLineColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->SetTitleColor(kBlack);
    axis->SetTextFont(42);
    axis->SetLabelFont(42);
    axis->CenterTitle();

    axis->Draw();


    


    // Here loop over the cloned histogram of nue, and zero out the flux after the threshold, this is so we can shade out the thresholded area
    for (unsigned int p=xbin_th; p < h_nue_clone->GetNbinsX()+1; p++){
        h_nue_clone->SetBinContent(p, 0);
    }

    h_nue_clone->SetLineColor(kBlue+2);
    h_nue_clone->SetLineWidth(0);
    h_nue_clone->SetFillColorAlpha(46, 0.4);
    if (!draw_averge)h_nue_clone->Draw("hist,same");
    // h_nuebar->Draw("hist,same");



    // Get the energy values
    double average_num = 0; // flux numerator
    double average_den = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1) +  h_nuebar->Integral(xbin_th, h_nuebar->GetNbinsX()+1); // flux denominator

    int nue_flux_bins  = h_nue->GetNbinsX();
    int anue_flux_bins = h_nuebar->GetNbinsX();

    for (int k = 1; k < nue_flux_bins+1; k++ ){
        if (k < xbin_th) continue; // skip the thresholded bins
        average_num = average_num + (h_nue->GetBinContent(k) * h_nue->GetBinCenter(k));
    }

    for (int k = 1; k < anue_flux_bins+1; k++ ){
        if (k < xbin_th) continue; // skip the thresholded bins
        average_num = average_num + (h_nuebar->GetBinContent(k) * h_nuebar->GetBinCenter(k));
    }

    double average = average_num / average_den;
    std::cout << "Average Energy: " << average << std::endl;
    

    std::cout << "Total Nue Flux: " << h_nue->Integral() << std::endl;
    std::cout << "Total Nuebar Flux: " << h_nuebar->Integral() << std::endl;
    std::cout << "Total Nue Flux with threshold: " << flux_int_thresh_nue << std::endl;
    std::cout << "Total Nuebar Flux with threshold: " << flux_int_thresh_nuebar << std::endl;

    double sum_flux_thresh_to_KDAR_end = h_nue->Integral( xbin_th, kdar_max) + h_nuebar->Integral( xbin_th, kdar_max);

    

    h_summed_flux_clone->SetLineWidth(0);
    h_summed_flux_clone->SetFillColorAlpha(46, 0.4);
    h_summed_flux_clone->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e} + #bar{#nu}_{e} / cm^{2} / 6 #times 10^{20} POT");
    if (draw_averge) h_summed_flux_clone->Draw("hist,same");

    
    double summed_flux_integral = summed_flux->Integral( xbin_th, summed_flux->GetNbinsX()+1);
    
    std::cout <<  "summed_flux_integral: " << summed_flux_integral << std::endl;
    std::cout << "Fraction of flux from thresh to KDAR end point: " << 100 * sum_flux_thresh_to_KDAR_end / summed_flux_integral << std::endl;

    // Get the bin where the average is located
    double average_bin = 0;
    average_bin  = summed_flux->GetXaxis()->FindBin(average);
    std::cout << "Average bin number: " << average_bin << std::endl;

    // Get the max error bar = 34% of the flux after the average
    double integral_max = 0;
    double max_bin_val = 0;
    for (int bin = average_bin; bin < summed_flux->GetNbinsX()+1; bin++){
        
        integral_max = summed_flux->Integral( average_bin, bin);
        
        if (integral_max / summed_flux_integral >= 0.34){
            max_bin_val = summed_flux->GetBinCenter(bin);
            break;
        }
    }
       
    // Get the min error bar = 34% of the flux beforethe average
    double min_sum = 0;
    double min_bin_val = 0;
    for (int bin = average_bin; bin > 0; bin--){
        if (bin < xbin_th) continue; // skip the thresholded bins
        
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
    double data_xsec = 0.667 * scale; // Scale is to take it to the flux axes

    double data_stat_err_plus = 0.144 * scale;
    double data_stat_err_minus = 0.144 * scale;

    double data_stat_plus_sys_err_plus = 0.267 * scale;
    double data_stat__plus_sys_err_minus = 0.267 * scale;

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
    TLegend *leg;
    if (!draw_averge)leg = new TLegend(0.17, 0.5, 0.48, 0.9);
    if (draw_averge)leg  = new TLegend(0.17, 0.5, 0.48, 0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (!draw_averge)leg->AddEntry(h_nue,             "NuMI #nu_{e} Flux",                       "f");
    if (!draw_averge)leg->AddEntry(h_nuebar,          "NuMI #bar{#nu}_{e} Flux",                 "f");
    if (!draw_averge)leg->AddEntry(h_spline_nue,      "GENIE #nu_{e} CC Cross-Section",          "l");
    if (!draw_averge)leg->AddEntry(h_spline_nuebar,   "GENIE #bar{#nu}_{e} CC Cross-Section",    "l");
    if (draw_averge)leg->AddEntry(summed_flux,        "NuMI #nu_{e} + #bar{#nu}_{e} Flux",       "f");
    if (draw_averge)leg->AddEntry(h_spline_average,   "GENIE #nu_{e} + #bar{#nu}_{e} CC Cross-Section",    "l");
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

    // summed_flux->GetXaxis()->SetTitle("test");
    // c->Modified();


    if (!draw_averge)c->Print("../../Analysis/plots/flux_combo_morebins_v3.pdf");
    if (draw_averge)c->Print("../../Analysis/plots/flux_combo_morebins_v3_average.pdf");


    // Now create the other flux integrated cross sec plot
    TCanvas *c2  = new TCanvas();
    double x1[]  = {0, 1};
    double y1[]  = {6.67e-39, 6.67e-39};
    double ex1[] = {0.0, 0.0};
    double ey1[] = {0.267e-38, 0.267e-38};
    double ey2[] = {0.144e-38, 0.144e-38};
    
    // Stats band
    auto g_xsec_sys = new TGraphErrors(2, x1, y1, ex1, ey1);
    g_xsec_sys->SetFillColor(kGray);
    g_xsec_sys->SetFillColorAlpha(12, 0.15);
    g_xsec_sys->SetTitle(";;#nu_{e} + #bar{#nu}_{e} CC Cross-Section [cm^{2}]");
    g_xsec_sys->GetXaxis()->SetRangeUser(0,1);
    g_xsec_sys->GetYaxis()->SetRangeUser(0.2e-38,1.4e-38);
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
    double y_xsec[]  = {6.67e-39, 6.67e-39};
    auto g_xsec = new TGraphErrors(2, x1, y_xsec);
    g_xsec->SetLineColor(kBlack);
    g_xsec->SetLineWidth(2);
    g_xsec->Draw("same");

    // Genie nue
    double y_nue[]  = {9.4e-39, 9.4e-39};
    auto g_xsec_nue = new TGraphErrors(2, x1, y_nue);
    g_xsec_nue->SetLineColor(kBlue+2);
    g_xsec_nue->SetLineWidth(2);
    g_xsec_nue->SetLineStyle(3);
    g_xsec_nue->Draw("same");

    // Genie nuebar
    double y_nuebar[]  = {2.85e-39, 2.85e-39};
    auto g_xsec_nuebar = new TGraphErrors(2, x1, y_nuebar);
    g_xsec_nuebar->SetLineColor(kGreen+2);
    g_xsec_nuebar->SetLineWidth(2);
    g_xsec_nuebar->SetLineStyle(3);
    g_xsec_nuebar->Draw("same");

    // Genie nue + nuebar
    double y_nue_nuebar[]  = {6.9e-39, 6.9e-39};
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
