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

    TH1D *h_spline_nue    = new TH1D("h_spline_nue", ";Energy [GeV]; Cross Section",num_bins, 0, 125);
    TH1D *h_spline_nuebar = new TH1D("h_spline_nuebar", ";Energy [GeV]; Cross Section",num_bins, 0, 125);

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
    h_nue->SetTitle(";Electron Neutrino Energy [GeV];#nu_{e}/#bar{#nu_{e}} / cm^{2} / 6e20 POT");
    h_nue->GetXaxis()->SetRangeUser(0,3);
    h_nue->GetYaxis()->SetRangeUser(0,150.0e6);
    h_nue->SetLineColor(kBlue+2);
    h_nue->SetFillColor(17);
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
    h_spline_nue->SetLineStyle(2);
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
    axis->SetTitle("#nu_{e}/#bar{#nu_{e}} CC Cross Section [1e-38 cm^{2}]");
    axis->SetTitleOffset(1.2);
    axis->SetLineColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->SetTitleColor(kBlack);
    
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
    gr->Draw("L,P,same");

    // Draw the Legend
    TLegend *leg = new TLegend(0.2, 0.6, 0.4, 0.8);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_nue,             "NuMI #nu_{e} Flux",                       "l");
    leg->AddEntry(h_nuebar,          "NuMI #bar{#nu_{e}} Flux",                 "l");
    leg->AddEntry(h_spline_nue,      "GENIE #nu_{e} CC Cross Section",          "l");
    leg->AddEntry(h_spline_nuebar,   "GENIE #bar{#nu_{e}} CC Cross Section",    "l");

    leg->AddEntry(gr, "Data #sigma_{#nu_{e}+#bar{#nu_{e}}} (stat+sys)",   "lep");

    leg->Draw();


}
