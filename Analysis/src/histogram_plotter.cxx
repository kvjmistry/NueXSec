#include "../include/histogram_plotter.h"

// -----------------------------------------------------------------------------
histogram_plotter::~histogram_plotter(){ 
    
    // Overwrite the destrctor of histogram helper
    std::cout << " " << std::endl;
}
// -----------------------------------------------------------------------------
void histogram_plotter::InitHistograms(){
}
// -----------------------------------------------------------------------------
void histogram_plotter::MakeStack(std::string hist_name, std::string cut_name, bool area_norm,  bool logy, const char* x_axis_name,
                                     double data_scale_factor, double y_scale_factor, double intime_scale_factor, double dirt_scale_factor, 
                                     const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, const char* print_name ){

    
    std::vector<TH1D*> hist(k_classifications_MAX);
    
    for (unsigned int i=0; i <classification_dirs.size(); i++){
        
        if (i == k_leg_data) 
            _utility_instance.GetHist(f_nuexsec, hist.at(i), Form("Data/Stack/%s/%s/%s_%s_%s", cut_name.c_str(), classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), classification_dirs.at(i).c_str()));
        if (i == k_leg_ext)
            _utility_instance.GetHist(f_nuexsec, hist.at(i), Form("EXT/Stack/%s/%s/%s_%s_%s",  cut_name.c_str(), classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), classification_dirs.at(i).c_str()));
        if (i == k_leg_dirt)
            _utility_instance.GetHist(f_nuexsec, hist.at(i), Form("Dirt/Stack/%s/%s/%s_%s_%s", cut_name.c_str(), classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), classification_dirs.at(i).c_str()));
        
        if (hist.at(i) != NULL && ( i == k_leg_data || i == k_leg_ext || i == k_leg_dirt)) continue;

        _utility_instance.GetHist(f_nuexsec, hist.at(i), Form("MC/Stack/%s/%s/%s_%s_%s", cut_name.c_str(), classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), classification_dirs.at(i).c_str()));
    
    }

    // If any of the histograms are null then return as this function will break
    for (unsigned int i=0; i < hist.size(); i++){
        if (hist.at(i) == NULL) return;
    }

    // Choose whether to use a pvalue
    const bool p_value = false;

    TCanvas* c       = new TCanvas();
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    
    topPad   ->SetBottomMargin(0.05);
    bottomPad->SetTopMargin(0.04);
    bottomPad->SetBottomMargin(0.25);
    bottomPad->SetGridy();
    topPad   ->Draw();
    bottomPad->Draw();
    topPad   ->cd();

    THStack * h_stack = new THStack();
    
    // Generic actions
    for (unsigned int i=0; i < hist.size(); i++){
        hist.at(i)->SetStats(kFALSE);
        
        if (i == k_leg_data) continue; // Dont scale the data 
        
        // Implement scaling here
        if (i == k_leg_ext){
            hist.at(i)->Scale(intime_scale_factor);
            continue;
        }

        if (i == k_leg_dirt){
            hist.at(i)->Scale(dirt_scale_factor);
            continue;
        }
        
        hist.at(i)->Scale(data_scale_factor);
    }
    
    // Customise the histogram
    hist.at(k_nue_cc)       ->SetFillColor(30);
    hist.at(k_nue_cc_mixed) ->SetFillColor(38);
    hist.at(k_numu_cc)      ->SetFillColor(28);
    hist.at(k_nc_pi0)       ->SetFillColor(36);
    hist.at(k_cosmic)       ->SetFillColor(1);
    hist.at(k_nc)           ->SetFillColor(46);
    hist.at(k_nue_cc_out_fv)->SetFillColor(kTeal);
    hist.at(k_nc_mixed)     ->SetFillColor(42);
    hist.at(k_unmatched)    ->SetFillColor(12);
   
    hist.at(k_leg_ext)      ->SetFillColor(41);
    hist.at(k_leg_ext)      ->SetFillStyle(3345);

    hist.at(k_leg_dirt)     ->SetFillColor(2);
    hist.at(k_leg_dirt)     ->SetFillStyle(3354);

    hist.at(k_leg_data)     ->SetMarkerStyle(20);
    hist.at(k_leg_data)     ->SetMarkerSize(0.5);


    double integral_data = hist.at(k_leg_data)->Integral();
    double integral_mc_ext{0.0};

    // Normalisation
    if (area_norm && integral_data != 0) {

        // Get the integral of the MC + EXT
        for (unsigned int i=0; i < hist.size(); i++){
            if (i == k_leg_data) continue; // Dont use the data
            integral_mc_ext += hist.at(i)->Integral();
        
        }
        
        // Check the integral of ON - EXT
        TH1D * h_data_scaling_clone = (TH1D*) hist.at(k_leg_data)->Clone("h_data_scaling_clone");
        h_data_scaling_clone->Add(  hist.at(k_leg_ext), -1);
        double integral_on_minus_off = h_data_scaling_clone->Integral();
        
        if (integral_on_minus_off == 0) {
            std::cout << "unable to area normalise" << std::endl;
            integral_on_minus_off = 1;
        }
        delete h_data_scaling_clone;

        for (unsigned int i=0; i < hist.size(); i++){
            if (i == k_leg_data) continue; // Dont use the data
            
            hist.at(i)->Scale(integral_on_minus_off / integral_mc_ext);
            hist.at(i)->Scale(1. / integral_data);
        
        }

    }

    // Add the histograms to the stack
    for (unsigned int i=0; i < hist.size(); i++){
        if (i == k_leg_data) continue; // Dont use the data
        h_stack->Add(hist.at(i));
    }

    const double y_maximum = std::max(hist.at(k_leg_data)->GetMaximum(), h_stack->GetMaximum());

    // Set the axis in the case of a log plot
    if (logy == true){
        TH1D * h_scale_axes = (TH1D*)hist.at(k_leg_data)->Clone("h_scale_axes");
        
        if(hist.at(k_nue_cc)->GetMinimum() != 0.0) {h_scale_axes->SetMinimum(hist.at(k_nue_cc)->GetMinimum() / 2.); }
        
        if(hist.at(k_nue_cc)->GetMinimum() == 0.0) {h_scale_axes->SetMinimum(hist.at(k_nue_cc)->GetMinimum() + 0.0001 / 2.); }
        
        h_scale_axes->SetMaximum(y_maximum * (y_scale_factor * 500));
        h_scale_axes->SetLineColor(0);
        h_scale_axes->SetFillColor(0);
        h_scale_axes->GetYaxis()->SetTitle("Entries [A.U.]");
        h_scale_axes->SetTitle(" ");
        h_scale_axes->GetXaxis()->SetTitle(" ");
        h_scale_axes->GetXaxis()->SetLabelSize(0);
        h_scale_axes->GetXaxis()->SetLabelFont(0); // Absolute font size in pixel (precision 3)
        h_scale_axes->Draw();
        h_stack->Draw("same hist");

        h_scale_axes->GetYaxis()->SetTitle("Entries [A.U.]");
    }

    // Set the axis in the case of a non-log plot
    if(!logy) {
        h_stack->SetMinimum(0);
        h_stack->SetMaximum(y_maximum * y_scale_factor);
        h_stack->Draw("hist");
    }
    
    // Set the y axis of the stack
    if(!area_norm) h_stack->GetYaxis()->SetTitle("Entries");
    else           h_stack->GetYaxis()->SetTitle("Entries [A.U.]");
    
    h_stack->GetYaxis()->SetTitleFont(45);
    h_stack->GetYaxis()->SetTitleSize(18);
    h_stack->GetYaxis()->SetTitleOffset(1.30);
    h_stack->GetXaxis()->SetLabelOffset(10);
    hist.at(k_leg_data)->Draw("same PE");

    // MC error histogram
    TH1D * h_error_hist = (TH1D*) hist.at(k_nue_cc)->Clone("h_error_hist");
    for (unsigned int i=0; i < hist.size(); i++){
        if (i == k_leg_data || i == k_nue_cc) continue; // Dont use the data
        h_error_hist->Add(hist.at(i), 1);
    }
    
    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2 hist same");

    // Set the legend
    TLegend *leg_stack = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
    //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    leg_stack->AddEntry(hist.at(k_nue_cc),          "#nu_{e} CC",        "f");
    leg_stack->AddEntry(hist.at(k_nue_cc_mixed),    "#nu_{e} CC Mixed",  "f");
    leg_stack->AddEntry(hist.at(k_nue_cc_out_fv),   "#nu_{e} CC OutFV",  "f");
    leg_stack->AddEntry(hist.at(k_cosmic),          "Cosmic",            "f");
    leg_stack->AddEntry(hist.at(k_numu_cc),         "#nu_{#mu} CC",      "f");
    leg_stack->AddEntry(hist.at(k_nc),              "NC",                "f");
    leg_stack->AddEntry(hist.at(k_nc_pi0),          "NC #pi^{0}",        "f");
    leg_stack->AddEntry(hist.at(k_nc_mixed),        "NC Mixed",          "f");
    leg_stack->AddEntry(hist.at(k_unmatched),       "Unmatched",         "f");
    leg_stack->AddEntry(hist.at(k_leg_dirt),        "Dirt",              "f");
    leg_stack->AddEntry(hist.at(k_leg_ext),         "InTime (EXT)",      "f");
    leg_stack->Draw();

    if(!logy) h_stack->GetYaxis()->SetRangeUser(0, y_maximum * y_scale_factor);
    else      topPad->SetLogy();
    
    // Calculate the chi2
    TH1D * h_last = (TH1D*) h_stack->GetStack()->Last();
    std::vector <double> chi2  = Chi2Calc(h_last, hist.at(k_leg_data), area_norm, integral_data);
    //chi2 : chi2/ndf, mc+ext, data

    
    // Plot the Reduced Chi2
    TPaveText * pt = new TPaveText(.46,.80,.73,1.06, "NBNDC");
    std::ostringstream o_string;
    o_string.precision(3);
    o_string << std::fixed;
    o_string << float(chi2.at(0));
    std::string convert_string = o_string.str();
    std::string chi2_string = "#chi_{Stat}^{2}/DOF=" + convert_string;
    pt->AddText(chi2_string.c_str());
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->Draw();

    // Num events
    TPaveText * pt2 = new TPaveText(.13,.80,.46,1.06, "NBNDC");
    std::ostringstream o_string2a;
    std::ostringstream o_string2b;
    if(!area_norm) {
        o_string2a << int(chi2.at(2));
        o_string2b << int(chi2.at(1));
    }
    if(area_norm) {
        o_string2a << int(chi2.at(2) * integral_data);
        o_string2b << int(chi2.at(1) * integral_data);
    }

    std::string convert_string2a = o_string2a.str();
    std::string convert_string2b = o_string2b.str();
    std::string chi2_string2 = "Data: " + convert_string2a + "|MC+EXT:" + convert_string2b;
    pt2->AddText(chi2_string2.c_str());
    pt2->SetFillStyle(0);
    pt2->SetBorderSize(0);
    // This is removed for public distributions
    pt2->Draw();

    // Num bins
    TPaveText * pt3 = new TPaveText(.60,.80,.73,.973, "NBNDC");
    std::ostringstream o_string3;
    o_string3 << int(chi2.at(3));
    std::string convert_string3 = o_string3.str();
    std::string ndf_string = "DOF=" + convert_string3;
    pt3->AddText(ndf_string.c_str());
    pt3->SetFillStyle(0);
    pt3->SetBorderSize(0);
    pt3->Draw();

    // p value -- optional
    TPaveText * pt4 = new TPaveText(.45,.80,.60,.973, "NBNDC");
    std::ostringstream o_string4;
    o_string4.precision(4);
    o_string4 << std::fixed;
    o_string4 << chi2.at(4);
    std::string convert_string4 = o_string4.str();
    std::string p_string = "P=" + convert_string4;
    pt4->AddText(p_string.c_str());
    pt4->SetFillStyle(0);
    pt4->SetBorderSize(0);
    if(p_value) {pt4->Draw(); }

    bottomPad->cd();

    // Now create the ratio of data to MC
    TH1D * ratioPlot = (TH1D*) hist.at(k_leg_data)->Clone("ratioPlot");
    TH1D * h_mc_ext_sum = (TH1D*) hist.at(k_nue_cc)->Clone("h_mc_ext_sum");
    
    for (unsigned int i=0; i < hist.size(); i++){
        if (i == k_leg_data || i == k_nue_cc ) continue; // Dont use the data and nue cc because already been cloned
        h_mc_ext_sum->Add(hist.at(i), 1);
    }
   
    ratioPlot->GetXaxis()->SetLabelSize(12);
    ratioPlot->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratioPlot->GetYaxis()->SetLabelSize(11);
    ratioPlot->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratioPlot->GetXaxis()->SetTitleOffset(3.0);
    ratioPlot->GetXaxis()->SetTitleSize(17);
    ratioPlot->GetXaxis()->SetTitleFont(46);
    ratioPlot->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
    ratioPlot->Add(h_mc_ext_sum, -1);
    ratioPlot->Divide(h_mc_ext_sum);
    ratioPlot->GetYaxis()->SetRangeUser(-1,1);
    ratioPlot->GetXaxis()->SetTitle(x_axis_name);
    ratioPlot->GetYaxis()->SetTitle("(Data - MC) / MC ");
    ratioPlot->GetYaxis()->SetTitleSize(13);
    ratioPlot->GetYaxis()->SetTitleFont(44);
    ratioPlot->GetYaxis()->SetTitleOffset(1.5);
    ratioPlot->SetTitle(" ");
    ratioPlot->Draw();

    // Now doing this stuff on the bottom pad
	//x_min, y_min, x_max, y_max
	// Reduced chi2
	TPaveText * pt_bottom = new TPaveText(.12, .80, .30, .96, "NBNDC");
	std::ostringstream o_string_bottom;
	o_string_bottom.precision(3);
	o_string_bottom << std::fixed;
	o_string_bottom << float(chi2.at(0) * chi2.at(3));
	std::string convert_string_bottom = o_string_bottom.str();

	std::ostringstream o_string3_bottom;
	o_string3_bottom << int(chi2.at(3));
	std::string convert_string3_bottom = o_string3_bottom.str();

	std::string chi2_string_bottom = "#chi_{Stat}^{2}/DOF=(" + convert_string_bottom + "/" + convert_string3_bottom + ")";
	pt_bottom->AddText(chi2_string_bottom.c_str());
	pt_bottom->SetFillStyle(0);
	pt_bottom->SetBorderSize(0);
	// pt_bottom->Draw();

    c->Print(print_name);


}
// -----------------------------------------------------------------------------
std::vector<double> histogram_plotter::Chi2Calc(TH1D * h_mc_ext, TH1D * h_data, const bool area_norm, const double return_norm){
    const int n_bins = h_mc_ext->GetNbinsX();

    const double f_1 = h_mc_ext->Integral();
    const double f_2 = h_data  ->Integral();

    //area normalised?
    TH1D * h_mc_ext_clone = (TH1D*)h_mc_ext->Clone("h_mc_ext_clone");
    TH1D * h_data_clone   = (TH1D*)h_data  ->Clone("h_data_clone");
    
    if(!area_norm) {h_mc_ext_clone->Scale(f_2/f_1); }
    
    if(area_norm) {
        //this keeps them area normalised,
        //but at the original values, not 0->1
        //which messes with the chi2 and p calc
        h_mc_ext_clone->Scale(return_norm);
        h_data_clone  ->Scale(return_norm);

        const double f_1_adj = h_mc_ext->Integral();
        const double f_2_adj = h_data->Integral();
        h_mc_ext_clone->Scale(f_2_adj/f_1_adj);

    }
    
    //h_data_clone->Scale(1./f_2);

    std::vector <double> chi2;
    double chi2_val     = 0;
    double n_mc_ext_val = 0;
    double n_data_val   = 0;
    
    // Loop over each bin
    for ( int i = 1; i < n_bins; i++) {
        const double n_mc_ext = h_mc_ext_clone->GetBinContent(i);
        const double n_data   = h_data_clone  ->GetBinContent(i);

        //don't calculate chi2 for bins where no comparison possible
        if(n_data == 0 || n_mc_ext == 0) { continue; }

        //chi2_val += (pow((n_mc_ext - n_data),2) / n_mc_ext);
        //chi2_val += (pow((n_data - n_mc_ext),2)) / (((n_data * f_2) / pow(f_2, 2)) + ((n_mc_ext * f_1) / pow(f_1, 2)));
        chi2_val += 2 * (n_mc_ext - n_data + (n_data * TMath::Log(n_data/n_mc_ext)));

        n_mc_ext_val += n_mc_ext;
        n_data_val   += n_data;
    }
    
    const double reduced_chi2 = chi2_val / (n_bins - 1);
    const double p = TMath::Prob(chi2_val, n_bins);

    chi2.push_back(reduced_chi2);
    //correct this value back to the un-normalised
    //chi2.push_back(n_mc_ext_val * f_1);
    //chi2.push_back(n_data_val * f_2);
    chi2.push_back(n_mc_ext_val * (f_1 / f_2));
    chi2.push_back(n_data_val);
    chi2.push_back(n_bins - 1);
    chi2.push_back(p);
    return chi2;
}