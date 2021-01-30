



void check_response(){

    // Open the file with the distributions
    TFile *f = TFile::Open("files/xsec_result_run1.root", "READ");
    // TFile *f = TFile::Open("../../Wiener-SVD-Unfolding/wiener_example1.root", "READ");
    

    // Load in the MC X sec as a function of true energy
    TH1D *h_mcxsec_true = (TH1D*)f->Get("elec_E/wiener/h_mc_xsec_true");
    // TH1D *h_mcxsec_true = (TH1D*)f->Get("htrue_signal");

    // response matrix
    TH2D *h_response = (TH2D*)f->Get("elec_E/wiener/h_response");
    // TH2D *h_response = (TH2D*)f->Get("hR");

    // MC event rate reco
    TH1D *h_mcer_reco = (TH1D*)f->Get("elec_E/er/h_mc_xsec_reco");
    // TH1D *h_mcer_reco = (TH1D*)f->Get("elec_E/wiener/h_data_xsec_stat_reco");

    // Clone histogram and clear it
    TH1D* h_mcxsec_smear = (TH1D*)h_mcer_reco->Clone();

    // h_mcxsec_true->Scale(1.0, "width");
    // h_mcer_reco->Scale(1.0, "width");
    
    // Clear the Bins
    for (int bin = 0; bin < h_mcxsec_smear->GetNbinsX()+2; bin++){
        h_mcxsec_smear->SetBinContent(bin, 0);
    }

    // for (int bin = 0; bin < h_mcxsec_true->GetNbinsX()+2; bin++){
    //     h_mcxsec_true->SetBinContent(bin, h_mcxsec_true->GetBinContent(bin)*h_mcxsec_true->GetBinWidth(bin));
    // }

    // Now do the multiplication
    for (int i=1; i < h_response->GetXaxis()->GetNbins()+1; i++){

        for (int j=1; j < h_response->GetYaxis()->GetNbins()+1; j++) { 
            
            h_mcxsec_smear->SetBinContent(i, h_mcxsec_smear->GetBinContent(i) + h_response->GetBinContent(i, j) * h_mcxsec_true->GetBinContent(j));
        }

    } 

    h_mcxsec_smear->Scale(1.0, "width");
    // h_mcer_reco->Scale(1.0, "width");

    TCanvas *c = new TCanvas();
    h_mcxsec_smear->SetLineWidth(2);
    h_mcxsec_smear->SetLineColor(kBlack);
    h_mcxsec_smear->Draw("hist");

    h_mcer_reco->SetLineColor(kRed);
    h_mcer_reco->SetLineWidth(1);
    h_mcer_reco->Draw("hist,same");





}