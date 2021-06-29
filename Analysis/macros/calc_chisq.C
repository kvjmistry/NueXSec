


void calc_chisq(){

    // Load in the TFile
    TFile* f = TFile::Open("files/nuexsec_mc_run1.root", "READ");
    TH1D* h = (TH1D*)f->Get("True/h_true_elec_gamma_MC_unselected");


    std::cout <<h->Integral()/h->GetNbinsX() << std::endl;

    double mean = 595.859;

    double chi = 0.0;
    for (int bin = 1; bin < h->GetNbinsX()+1; bin++){
        chi+= (h->GetBinContent(bin) - mean) * (1.0/h->GetBinError(bin) ) *(h->GetBinContent(bin) - mean);
    }

    std::cout << "Chsq/N: " <<chi << "/" << h->GetNbinsX() << ": " << chi/double(h->GetNbinsX())<< std::endl;


}
