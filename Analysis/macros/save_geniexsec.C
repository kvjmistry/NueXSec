


void save_geniexsec(){


    // Load in the Cross section prediction
    TFile *f_in = TFile::Open("../files/crosssec_run1.root", "READ");

    TH1D* h_E    = (TH1D*)f_in->Get("geniev3/true_el_E/h_run1_CV_0_true_el_E_mc_xsec");
    TH1D* h_cang = (TH1D*)f_in->Get("geniev3/true_el_cang/h_run1_CV_0_true_el_cang_mc_xsec");

    h_E->SetDirectory(0);
    h_cang->SetDirectory(0);

    f_in->Close();

    h_E->Scale(1.0, "width");
    h_cang->Scale(1.0, "width");

    TFile* f_out = TFile::Open("../files/xsec_result_run1_paper_forandy.root", "UPDATE");

    h_E->Write("mc_xsec_true_genie_v3_0_6_energy_new");
    h_cang->Write("mc_xsec_true_genie_v3_0_6_angle_new");


    f_out->Close();
}