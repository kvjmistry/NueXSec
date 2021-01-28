

// Script to make the efficiency plot in MCC8. Uses hard coded values unfortunantely

void make_efficiency_plot(){


    TFile *f_mc;

    //std::vector<double> efficiency_v = {1, 0.5489, 0.38, 0.315, 0.261, 0.132, 0.0903}; // efficiency vector -- with mixed as bkg
    //std::vector<double> purity_v = {0, 0.00736, 0.03569, 0.0418, 0.1514, 0.2928, 0.3973}    ; // purity vector -- with mixec as bkg
    //std::vector<double> purity_nu_v = {0, 0.1225, 0.1612, 0.1943, 0.3253, 0.5142, 0.6978}    ; // purity nu Only vector -- with mixed as bkg
    std::vector<double> efficiency_v = {1, 0.694, 0.451, 0.360, 0.299, 0.15, 0.091}; // efficiency vector -- with mixed as sig
    std::vector<double> purity_v = {0, 0.0092, 0.042, 0.047, 0.172, 0.331, 0.385}    ; // purity vector -- with mixec as sig
    std::vector<double> purity_nu_v = {0, 0.077, 0.162,  0.184, 0.352, 0.569, 0.669}    ; // purity nu Only vector -- with mixed as sig
    std::vector<std::string> names = {"No Selection (0)","Pre-selection (1)", "Flash Matching (2)", "Vertex Reco. Quality (3)", "Shower Hit Threshold (4)", "Electron-like Shower (5)", "Final Tuning (6)"};

    double efficiency, purity;

    TCanvas *c = new TCanvas();
    TH1D* h_eff = new TH1D("h_efficiency", "", efficiency_v.size(), 0, efficiency_v.size());
    TH1D* h_pur = new TH1D("h_purity", "", efficiency_v.size(), 0, efficiency_v.size());
    TH1D* h_pur_nu = new TH1D("h_purity_nu", "", efficiency_v.size(), 0, efficiency_v.size());

    // c->SetGrid();
    c->SetGridy();

    TLegend *leg_stack = new TLegend(0.54, 0.89, 0.84, 0.70);
    leg_stack->SetBorderSize(0);
    //leg_stack->SetFillStyle(0);


    for (unsigned int k=0; k < efficiency_v.size();k++){
        h_eff ->Fill(names.at(k).c_str(), efficiency_v.at(k));
        h_pur ->Fill(names.at(k).c_str(), purity_v.at(k));
        h_pur_nu ->Fill(names.at(k).c_str(), purity_nu_v.at(k));
        h_eff->SetBinError(k+1, 0);
        h_pur->SetBinError(k+1, 0);
        h_pur_nu->SetBinError(k+1, 0);
    }
    
    leg_stack->AddEntry(h_eff, "Efficiency","lp");
    leg_stack->AddEntry(h_pur, "Purity",    "lp");
    leg_stack->AddEntry(h_pur_nu, "Purity (Beam Only)",    "lp");

    h_eff->GetYaxis()->SetLabelSize(0.05);
    h_eff->GetYaxis()->SetRangeUser(0,1.1);
    h_eff->GetXaxis()->SetLabelSize(0.05);
    h_eff->GetXaxis()->SetLabelOffset(0.01);
    h_eff->SetStats(kFALSE);
    h_eff->SetMarkerStyle(20);
    h_eff->SetMarkerSize(0.5);
    h_eff->SetLineWidth(2);

    h_eff->Draw("LP");

    h_pur->SetLineColor(kRed+2);
    h_pur->SetStats(kFALSE);
    h_pur->SetMarkerStyle(20);
    h_pur->SetMarkerSize(0.5);
    h_pur->SetLineWidth(2);
    h_pur->Draw("LP,same");

    h_pur_nu->SetLineColor(kGreen+2);
    h_pur_nu->SetStats(kFALSE);
    h_pur_nu->SetMarkerStyle(20);
    h_pur_nu->SetMarkerSize(0.5);
    h_pur_nu->SetLineWidth(2);
    h_pur_nu->Draw("LP,same");


    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.15);

    // Draw vertical lines to help the eye
    TLine *line;
    for (unsigned int l=1; l < efficiency_v.size()+1; l++){
        line  = new TLine( h_eff->GetBinCenter(l) ,   0 , h_eff->GetBinCenter(l)  ,  1.1);
        line->SetLineColor(12);
        line->SetLineStyle(kDotted);
        line->Draw();
    }

    h_eff->GetXaxis()->SetTickLength(0.00);
    h_pur->GetXaxis()->SetTickLength(0.00);

    TPaveText *pt = new TPaveText(0.3, 0.875, 0.3, 0.875,"NDC");
    pt->AddText("MicroBooNE Simulation");
    pt->SetTextColor(kBlack);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.04);
    pt->Draw();

    leg_stack->Draw();



    leg_stack->Draw();
   
    //TPaveText *pt;
    //pt = new TPaveText(0.21, 0.88, 0.21, 0.88,"NDC");
    //pt->AddText("MicroBooNE");
    //pt->SetTextColor(kBlack);
    //pt->SetBorderSize(0);
    //pt->SetFillColor(0);
    //pt->SetFillStyle(0);
    //pt->SetTextSize(0.04);
    //pt->Draw();


    c->Print("efficiency_plot_mcc8.pdf");




}
