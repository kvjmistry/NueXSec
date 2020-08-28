// enum to switch file type 
enum type {k_mc, k_data, k_ext, k_dirt, k_type_MAX}; 

void Draw_Data_MC_Ratio(TCanvas *c, double ratio, double x1, double y1, double x2, double y2){
    c->cd();

    // 0.34, 0.936, 0.34, 0.936

    TPaveText *pt;

    pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->AddText(Form("Data/MC Ratio: %2.2f", ratio));
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}

void MakeFlashPlot(const char *print_name, std::string histname, std::vector<TFile*> file_v) {

    std::vector<TH1D *> hist(k_type_MAX);
    std::vector<double> hist_integrals(k_type_MAX, 0.0); // The integrals of all the histograms
    double integral_mc_ext = 0.0;


    bool area_norm = false;

    TH1D *h_ratio;
    TH1D *h_ratio_error;
    TH1D *h_mc_ext_sum;

    TPad *topPad;
    TPad *bottomPad;
    TCanvas * c = new TCanvas(Form("c_%s", histname.c_str()), "c", 500, 500);
    THStack *h_stack = new THStack();

    for (int k = 0; k < k_type_MAX; k++)
    {
        file_v.at(k)->cd();
        hist.at(k) =  (TH1D*)file_v.at(k)->Get(histname.c_str());
        if (hist.at(k) == NULL)
        {
            std::cout << "Couldn't get all the flash histograms so exiting function..." << std::endl;
            return;
        }
    }

    topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    topPad->SetBottomMargin(0.05);
    topPad->SetTopMargin(0.15);
    bottomPad->SetTopMargin(0.04);
    bottomPad->SetBottomMargin(0.25);
    bottomPad->SetGridy();
    topPad->SetLeftMargin(0.15);
    topPad->SetRightMargin(0.1);
    bottomPad->SetLeftMargin(0.15);
    bottomPad->SetRightMargin(0.1);
    topPad->Draw();
    bottomPad->Draw();
    topPad->cd();

    for (unsigned int i = 0; i < hist.size(); i++)
    {

        if (i == k_data)
        {

            hist.at(i)->SetStats(kFALSE);
            hist.at(k_data)->SetMarkerStyle(20);
            hist.at(k_data)->SetMarkerSize(0.5);
            hist_integrals.at(k_data) = hist.at(k_data)->Integral();
        }

        // Scale EXT
        else if (i == k_ext)
        {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(1.0154*0.98); 
            // hist.at(i)->Scale(1.0154);
            hist.at(k_ext)->SetFillColor(41);
            hist.at(k_ext)->SetFillStyle(3345);
            hist_integrals.at(k_ext) = hist.at(k_ext)->Integral();
            integral_mc_ext += hist_integrals.at(k_ext);
        }

        // Scale Dirt
        else if (i == k_dirt)
        {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(0.16411);
            hist.at(k_dirt)->SetFillColor(2);
            hist.at(k_dirt)->SetFillStyle(3354);
            hist_integrals.at(k_dirt) = hist.at(k_dirt)->Integral();
            integral_mc_ext += hist_integrals.at(k_dirt);
        }

        // Scale MC
        else
        {
            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(0.1301);
            hist.at(i)->SetFillColor(30);
            hist_integrals.at(k_mc) = hist.at(k_mc)->Integral();
            integral_mc_ext += hist_integrals.at(k_mc);
        }
    }

    // Normalisation by area
    if (area_norm)
    {

        if (integral_mc_ext != 0)
        {

            for (unsigned int i = 0; i < hist.size(); i++)
            {
                if (i == k_data)
                    continue; // Dont scale the data
                // if (i == 0) std::cout << "area norm scale factor: "  << hist_integrals.at(k_plot_data) / integral_mc_ext << std::endl;
                hist.at(i)->Scale(hist_integrals.at(k_data) / integral_mc_ext);
            }
        }
    }

    // Add the histograms to the stack
    h_stack->Add(hist.at(k_ext));
    h_stack->Add(hist.at(k_mc));
    h_stack->Add(hist.at(k_dirt));

    h_stack->Draw("hist");
    hist.at(k_data)->Draw("same PE");

    if (area_norm)
        h_stack->GetYaxis()->SetTitle("Entries A.U. ");
    else
        h_stack->GetYaxis()->SetTitle("Entries");

    h_stack->GetYaxis()->SetTitleSize(0.05);
    h_stack->GetYaxis()->SetLabelSize(0.05);
    h_stack->GetXaxis()->SetLabelSize(0);
    h_stack->GetXaxis()->SetRangeUser(2, 23);
    if (histname == "h_flash_singlebin") h_stack->GetXaxis()->SetRangeUser(5.6,15.4);

    // MC error histogram ------------------------------------------------------
    TH1D *h_error_hist = (TH1D *)hist.at(k_mc)->Clone("h_error_hist");

    for (unsigned int i = 0; i < hist.size(); i++)
    {
        if (i == k_data)
            continue; // Dont use the data
        if (i == k_mc)
            continue; // Aleady got this histogram from the clone

        h_error_hist->Add(hist.at(i), 1);
    }

    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2, same");

    TLegend *leg_stack;
    if (histname == "h_flash_singlebin") leg_stack = new TLegend(0.7, 0.2, 0.9, 0.45);
    else leg_stack = new TLegend(0.7, 0.6, 0.9, 0.85);
    leg_stack->SetBorderSize(0);
    if (histname != "h_flash_singlebin")leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(hist.at(k_data), "Beam-On Data", "lep");
    leg_stack->AddEntry(hist.at(k_dirt), "Out-of Cryo", "f");
    leg_stack->AddEntry(hist.at(k_mc),   "MC", "f");
    leg_stack->AddEntry(hist.at(k_ext),  "Beam-Off Data", "f");

    leg_stack->Draw();

    bottomPad->cd();

    h_ratio = (TH1D *)hist.at(k_data)->Clone("h_ratio");
    h_mc_ext_sum = (TH1D *)hist.at(k_mc)->Clone("h_mc_ext_sum");

    for (unsigned int i = 0; i < hist.size(); i++)
    {
        if (i == k_data || i == k_mc)
            continue; // Dont use the data and nue cc because already been cloned
        h_mc_ext_sum->Add(hist.at(i), 1);
    }

    // h_ratio->Add(h_mc_ext_sum, -1);
    h_ratio->Divide(h_mc_ext_sum);

    h_ratio->GetXaxis()->SetLabelSize(0.13);
    h_ratio->GetXaxis()->SetTitleOffset(0.9);
    h_ratio->GetXaxis()->SetTitleSize(0.13);
    h_ratio->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

    // For percent difference
    // h_ratio->GetYaxis()->SetTitle("(Data - MC) / MC ");
    // h_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);

    // For ratio

    h_ratio->GetYaxis()->SetRangeUser(0.80, 1.20);
    h_ratio->GetXaxis()->SetRangeUser(2, 23);
    if (histname == "h_flash_singlebin") h_ratio->GetXaxis()->SetRangeUser(5.6,15.4);
    h_ratio->GetYaxis()->SetTitle("#frac{Beam-On}{(MC + Beam-Off)}");

    h_ratio->GetYaxis()->SetLabelSize(0.13);
    h_ratio->GetYaxis()->SetTitleSize(0.07);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->SetTitle(" ");
    h_ratio->GetYaxis()->CenterTitle();
    h_ratio->Draw("E");

    // Draw the error hist
    h_ratio_error = (TH1D *)h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");

    // Draw the run period on the plot
    // _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915, run_period);

    // Draw Data to MC ratio
    Draw_Data_MC_Ratio(c, double(hist_integrals.at(k_data) * 1.0 / integral_mc_ext * 1.0), 0.34, 0.936, 0.34, 0.936);

    // Draw other data specifc quantities
    // _util.Draw_Data_POT(c, Data_POT, 0.45, 0.915, 0.45, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);

    // if (area_norm)
    //     Draw_Area_Norm(c);

    c->Print(print_name);
}



void plot_flash(){

    std::vector<TFile*> file_v;
    file_v.resize(k_type_MAX);

    file_v.at(k_mc)   = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_mc_mcc8.root");
    file_v.at(k_data) = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_data_mcc8.root");
    file_v.at(k_ext)  = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_ext_mcc8.root");
    file_v.at(k_dirt) = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_dirt_mcc8.root");

    MakeFlashPlot("plots/h_flash.pdf", "h_flash", file_v);
    MakeFlashPlot("plots/h_flash_singlebin.pdf", "h_flash_singlebin", file_v);


}