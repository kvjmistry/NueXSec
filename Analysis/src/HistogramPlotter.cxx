#include "../include/HistogramPlotter.h"

// -----------------------------------------------------------------------------
HistogramPlotter::~HistogramPlotter(){

    // Make sure the file is closed
    // f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeHistograms(Utility _utility) {

    std::cout << "Creating histograms and making plots" << std::endl;

    _util = _utility;

    double Data_POT;

    area_norm = _util.area_norm;

    // Set the scale factors
    if (strcmp(_util.run_period, "1") == 0){
        Data_POT = _util.config_v.at(_util.k_Run1_Data_POT); // Define this variable here for easier reading
    }
    else if (strcmp(_util.run_period, "3") == 0){
        Data_POT = _util.config_v.at(_util.k_Run3_Data_POT); // Define this variable here for easier reading
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
        exit(1);
    }

    // Set the variation mode
    if (std::string(_util.variation) != "empty") std::cout << "Using varmode for detector variation:\t " << _util.variation << std::endl;

    Initalise();

    // Only do this stuff for the CV -- unless we really want these plots for each detvar, then we need to change the file paths!!
    if (std::string(_util.variation) == "empty") {

        // Create a set of strings for creating a dynamic directory
        // Directory structure that is created will take the form plots/<cut>/
        _util.CreateDirectory("Flash");

        // Flash time plots
        MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_time.pdf", _util.run_period), "h_flash_time");
        MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_time_single_bin.pdf", _util.run_period), "h_flash_time_single_bin"); // Single bin in the beam window
        MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_time_sid1.pdf", _util.run_period), "h_flash_time_sid1"); // Slice ID 1
        MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_time_sid0.pdf", _util.run_period), "h_flash_time_sid0"); // Slice ID 0

        // Flash PE plots
        MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_pe.pdf", _util.run_period), "h_flash_pe");
        MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_pe_sid1.pdf", _util.run_period), "h_flash_pe_sid1"); // Slice ID 1
        MakeFlashPlot(Data_POT, Form("plots/run%s/Flash/flash_pe_sid0.pdf", _util.run_period), "h_flash_pe_sid0"); // Slice ID 0

        // On beam minus off beam plots

        // Flash time
        MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_time_OMO.pdf", _util.run_period), "h_flash_time");
        MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_time_OMO_sid1.pdf", _util.run_period), "h_flash_time_sid1"); // Slice ID 1
        MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_time_OMO_sid0.pdf", _util.run_period), "h_flash_time_sid0"); // Slice ID 0

        // Flash PE plots
        MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_pe_OMO.pdf", _util.run_period), "h_flash_pe");
        MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_pe_OMO_sid1.pdf", _util.run_period), "h_flash_pe_sid1"); // Slice ID 1
        MakeFlashPlotOMO(Data_POT, Form("plots/run%s/Flash/flash_pe_OMO_sid0.pdf", _util.run_period), "h_flash_pe_sid0"); // Slice ID 0

        // Create the Efficiency Folder
        _util.CreateDirectory("Efficiency");

        MakeEfficiencyPlot(Form("plots/run%s/Efficiency/Integrated_Efficiency_Purity.pdf", _util.run_period));

        MakeEfficiencyPlotByCut("h_true_nu_E", false);
        MakeEfficiencyPlotByCut("h_true_elec_E", false);
        MakeEfficiencyPlotByCut("h_true_elec_E_rebin", true);
        MakeEfficiencyPlotByCut("h_true_nu_E_single_bin", true);
        MakeEfficiencyPlotByCut("h_true_nu_E_nue", false);
        MakeEfficiencyPlotByCut("h_true_nu_E_nuebar", false);
        MakeEfficiencyPlotByCut("h_true_nu_E_nue_single_bin", true);
        MakeEfficiencyPlotByCut("h_true_nu_E_nuebar_single_bin", true);

        // Create the interaction folder
        _util.CreateDirectory("Interaction");

        // Interaction Plot
        MakeInteractionPlot(Form("plots/run%s/Interaction/True_nue_nuebar_e_interaction_unselected.pdf", _util.run_period), "unselected","nue_nuebar", true);
        MakeInteractionPlot(Form("plots/run%s/Interaction/True_nue_nuebar_e_interaction_selected.pdf", _util.run_period), "selected","nue_nuebar", true);
        MakeInteractionPlot(Form("plots/run%s/Interaction/True_nue_e_interaction_unselected.pdf", _util.run_period), "unselected","nue", true);
        MakeInteractionPlot(Form("plots/run%s/Interaction/True_nue_e_interaction_selected.pdf", _util.run_period), "selected","nue", true);
        MakeInteractionPlot(Form("plots/run%s/Interaction/True_nuebar_e_interaction_unselected.pdf", _util.run_period), "unselected","nuebar", true);
        MakeInteractionPlot(Form("plots/run%s/Interaction/True_nuebar_e_interaction_selected.pdf", _util.run_period), "selected","nuebar", true);
        MakeInteractionEfficiency(Form("plots/run%s/Interaction/True_nue_e_interaction_efficiency.pdf", _util.run_period));

        // Create the 2D folder
        _util.CreateDirectory("2D");

        Plot2D_Signal_Background(Form("plots/run%s/2D/reco_shr_dEdx_shr_dist.pdf", _util.run_period), "h_reco_shr_dEdx_shr_dist");
        Plot2D_Signal_Background(Form("plots/run%s/2D/reco_shr_dEdx_shr_dist_post.pdf", _util.run_period), "h_reco_shr_dEdx_shr_dist_post");
        Plot2D_Signal_Background(Form("plots/run%s/2D/reco_shr_dEdx_max_shr_dist.pdf", _util.run_period), "h_reco_shr_dEdx_max_shr_dist");
        Plot2D_Signal_Background(Form("plots/run%s/2D/reco_shr_dEdx_max_shr_dist_post.pdf", _util.run_period), "h_reco_shr_dEdx_max_shr_dist_post");
        Plot2D_Signal_Background(Form("plots/run%s/2D/reco_shr_dEdx_shr_dist_large_dedx.pdf", _util.run_period), "h_reco_shr_dEdx_shr_dist_large_dedx");

        Plot2D_Signal_Background(Form("plots/run%s/2D/reco_shr_dEdx_moliere.pdf", _util.run_period), "h_reco_shr_dEdx_moliere");
        Plot2D_Signal_Background(Form("plots/run%s/2D/reco_shr_moliere_shr_dist.pdf", _util.run_period), "h_reco_shr_moliere_shr_dist");

        // Create the Truth folder
        _util.CreateDirectory("Truth");

        // Make 1D Histtograms
        // Save1DHists(Form("plots/run%s/Truth/true_nue_theta.pdf", _util.run_period), "h_true_nue_theta", _util.run_period);
        for (int i = 0; i < 2; i++){
            std::string cut_type = "unselected";
            if (i == 1) cut_type = "selected";

            Save1DHists(Form("plots/run%s/Truth/true_nue_theta_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_nue_theta", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_phi_%s.pdf",   _util.run_period, cut_type.c_str()), "h_true_nue_phi", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_theta_numi_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_nue_theta_numi", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_phi_numi_%s.pdf",   _util.run_period, cut_type.c_str()), "h_true_nue_phi_numi", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_angle_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_nue_angle", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_px_%s.pdf",    _util.run_period, cut_type.c_str()), "h_true_nue_px", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_py_%s.pdf",    _util.run_period, cut_type.c_str()), "h_true_nue_py", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_pz_%s.pdf",    _util.run_period, cut_type.c_str()), "h_true_nue_pz", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_e_%s.pdf",     _util.run_period, cut_type.c_str()), "h_true_nue_e", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nue_p_%s.pdf",     _util.run_period, cut_type.c_str()), "h_true_nue_p", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_vtx_x_%s.pdf",     _util.run_period, cut_type.c_str()), "h_true_vtx_x", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_vtx_y_%s.pdf",     _util.run_period, cut_type.c_str()), "h_true_vtx_y", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_vtx_z_%s.pdf",     _util.run_period, cut_type.c_str()), "h_true_vtx_z", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_vtx_x_sce_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_vtx_x_sce", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_vtx_y_sce_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_vtx_y_sce", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_vtx_z_sce_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_vtx_z_sce", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_nu_ang_targ_%s.pdf",   _util.run_period, cut_type.c_str()), "h_true_nu_ang_targ",   cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_elec_ang_targ_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_elec_ang_targ", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_elec_E_%s.pdf",     _util.run_period, cut_type.c_str()),     "h_true_elec_E", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_elec_phi_%s.pdf",   _util.run_period, cut_type.c_str()),     "h_true_elec_phi", cut_type, true);
            Save1DHists(Form("plots/run%s/Truth/true_elec_theta_%s.pdf", _util.run_period, cut_type.c_str()),     "h_true_elec_theta", cut_type, true);

            // Make the 2D histograms
            Save2DHists(Form("plots/run%s/Truth/h_true_nue_phi_theta_%s.pdf",          _util.run_period, cut_type.c_str()), "h_true_nue_phi_theta", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_nue_energy_theta_%s.pdf",       _util.run_period, cut_type.c_str()), "h_true_nue_energy_theta", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_nue_energy_phi_%s.pdf",         _util.run_period, cut_type.c_str()), "h_true_nue_energy_phi", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_nue_energy_angle_%s.pdf",       _util.run_period, cut_type.c_str()), "h_true_nue_energy_angle", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_nue_vtx_z_y_%s.pdf",            _util.run_period, cut_type.c_str()), "h_true_nue_vtx_z_y", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_nue_vtx_z_y_sce_%s.pdf",        _util.run_period, cut_type.c_str()), "h_true_nue_vtx_z_y_sce", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_elec_E_reco_elec_E_%s.pdf",     _util.run_period, cut_type.c_str()), "h_true_elec_E_reco_elec_E", cut_type, true);
            Save2DHists(Form("plots/run%s/Truth/h_true_nu_E_reco_nu_E_%s.pdf",         _util.run_period, cut_type.c_str()), "h_true_nu_E_reco_nu_E", cut_type, true);
            Save2DHists(Form("plots/run%s/Truth/h_true_elec_E_reco_elec_E_extra_bins_%s.pdf",     _util.run_period, cut_type.c_str()), "h_true_elec_E_reco_elec_E_extra_bins", cut_type, true);
            Save2DHists(Form("plots/run%s/Truth/h_true_nu_E_reco_nu_E_extra_bins_%s.pdf",         _util.run_period, cut_type.c_str()), "h_true_nu_E_reco_nu_E_extra_bins", cut_type, true);
            Save2DHists(Form("plots/run%s/Truth/h_true_nu_vtx_x_reco_nu_vtx_x_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_nu_vtx_x_reco_nu_vtx_x", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_nu_vtx_y_reco_nu_vtx_y_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_nu_vtx_y_reco_nu_vtx_y", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_nu_vtx_z_reco_nu_vtx_z_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_nu_vtx_z_reco_nu_vtx_z", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_shr_energy_purity_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_shr_energy_purity", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_shr_energy_completeness_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_shr_energy_completeness", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_shr_energy_resolution_reco_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_shr_energy_resolution_reco", cut_type, false);
            Save2DHists(Form("plots/run%s/Truth/h_true_shr_energy_resolution_true_%s.pdf", _util.run_period, cut_type.c_str()), "h_true_shr_energy_resolution_true", cut_type, false);

            // Normalised by reco (row)
            Save2DHistsNorm(Form("plots/run%s/Truth/h_true_elec_E_reco_elec_E_%s_row_norm_reco.pdf",     _util.run_period, cut_type.c_str()), "h_true_elec_E_reco_elec_E", cut_type, true, "reco");

            // Normalised by true (col)
            Save2DHistsNorm(Form("plots/run%s/Truth/h_true_elec_E_reco_elec_E_%s_col_norm_true.pdf",     _util.run_period, cut_type.c_str()), "h_true_elec_E_reco_elec_E", cut_type, true, "true");
        }

        // Stacked histograms for pi0
        // Create the Truth folder
        _util.CreateDirectory("pi0");

        // pi0 mass peak unweighted
        MakeStack("h_pi0_mass", " ",
                area_norm, false, 1.0, "Mass [MeV]",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/pi0/pi0_mass.pdf", _util.run_period), false, "classifications_pi0", false, false, true);

        // pi0 mass normlaisation fixed
        MakeStack("h_pi0_mass_norm", " ",
                area_norm, false, 1.0, "Mass [MeV] (Normalisation Fixed)",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/pi0/pi0_mass_norm.pdf", _util.run_period), false, "classifications_pi0", false, false, true);

        // pi0 mass peak unweighted
        MakeStack("h_pi0_mass_EScale", " ",
                area_norm, false, 1.0, "Mass [MeV] (E Dependent Scaling)",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/pi0/pi0_mass_EScale.pdf", _util.run_period), false, "classifications_pi0", false, false, true);

        // pi0 energy unweighted
        MakeStack("h_pi0_energy", " ",
                area_norm, false, 1.5, "#pi^{0} Energy [MeV]",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/pi0/pi0_energy.pdf", _util.run_period), false, "classifications_pi0", false, false, true);

        // pi0 energy normlaisation fixed
        MakeStack("h_pi0_energy_norm", " ",
                area_norm, false, 1.5, "#pi^{0} Energy [MeV] (Normalisation Fixed)",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/pi0/pi0_energy_norm.pdf", _util.run_period), false, "classifications_pi0", false, false, true);

        // pi0 energy unweighted
        MakeStack("h_pi0_energy_EScale", " ",
                area_norm, false, 1.5, "#pi^{0} Energy [MeV] (E Dependent Scaling)",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/pi0/pi0_energy_EScale.pdf", _util.run_period), false, "classifications_pi0", false, false, true);

        // Stacked histograms for numu
        // Create the Truth folder
        _util.CreateDirectory("numu");

        // Track Theta
        MakeStack("h_track_theta", " ",
                area_norm, false, 1.0, "Leading Track Theta [deg]",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/numu/track_theta.pdf", _util.run_period), false, "classifications_numu", false, false, true);

        // Track phi
        MakeStack("h_track_phi", " ",
                area_norm, false, 1.0, "Leading Track Phi [deg]",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/numu/track_phi.pdf", _util.run_period), false, "classifications_numu", false, false, true);

        // Track Cos theta
        MakeStack("h_track_cos_theta", " ",
                area_norm, false, 1.0, "Leading Track Cos(#theta)",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/numu/track_cos_theta.pdf", _util.run_period), false, "classifications_numu", false, false, true);

        // Track Cos theta
        MakeStack("h_muon_topo_score", " ",
                area_norm, false, 1.0, "Topological Score",  0.35, 0.85, 0.55, 0.85, Data_POT,
                Form("plots/run%s/numu/topo_score.pdf", _util.run_period), false, "classifications_numu", false, false, true);

    }

    // Loop over the cuts and plot histograms by plot type
    for (unsigned int i = 0; i < _util.k_cuts_MAX; i++) {

        // Create the directories
        if (std::string(_util.variation) == "empty")
            _util.CreateDirectory("/cuts/" + _util.cut_dirs.at(i));
        else
            _util.CreateDirectory("/detvar/" + std::string(_util.variation) + "/cuts/" + _util.cut_dirs.at(i));

        // Call the Make stack function for all the plots we want
        CallMakeStack( i, Data_POT);
    }

    

}
// -----------------------------------------------------------------------------
void HistogramPlotter::Initalise() {

    std::cout << "Initalising Histogram Plotter..." << std::endl;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject(_util.hist_file_name)) {
        f_nuexsec = TFile::Open(_util.hist_file_name);
    }
    else {
        std::cout << "Can't find histogram file!! " << __PRETTY_FUNCTION__ << std::endl;
        exit(1);
    }

}
// -----------------------------------------------------------------------------
std::vector<double> HistogramPlotter::Chi2Calc(TH1D *h_mc_ext, TH1D *h_data, const bool area_norm, const double return_norm) {
    const int n_bins = h_mc_ext->GetNbinsX();

    const double f_1 = h_mc_ext->Integral();
    const double f_2 = h_data->Integral();

    //area normalised?
    TH1D *h_mc_ext_clone = (TH1D *)h_mc_ext->Clone("h_mc_ext_clone");
    TH1D *h_data_clone = (TH1D *)h_data->Clone("h_data_clone");

    if (!area_norm) {
        h_mc_ext_clone->Scale(f_2 / f_1);
    }

    if (area_norm) {
        //this keeps them area normalised,
        //but at the original values, not 0->1
        //which messes with the chi2 and p calc
        h_mc_ext_clone->Scale(return_norm);
        h_data_clone->Scale(return_norm);

        const double f_1_adj = h_mc_ext->Integral();
        const double f_2_adj = h_data->Integral();
        h_mc_ext_clone->Scale(f_2_adj / f_1_adj);
    }

    //h_data_clone->Scale(1./f_2);

    std::vector<double> chi2;
    double chi2_val = 0;
    double n_mc_ext_val = 0;
    double n_data_val = 0;

    // Loop over each bin
    for (int i = 1; i < n_bins; i++) {
        const double n_mc_ext = h_mc_ext_clone->GetBinContent(i);
        const double n_data = h_data_clone->GetBinContent(i);

        //don't calculate chi2 for bins where no comparison possible
        if (n_data == 0 || n_mc_ext == 0) {
            continue;
        }

        //chi2_val += (pow((n_mc_ext - n_data),2) / n_mc_ext);
        //chi2_val += (pow((n_data - n_mc_ext),2)) / (((n_data * f_2) / pow(f_2, 2)) + ((n_mc_ext * f_1) / pow(f_1, 2)));
        chi2_val += 2 * (n_mc_ext - n_data + (n_data * TMath::Log(n_data / n_mc_ext)));

        n_mc_ext_val += n_mc_ext;
        n_data_val += n_data;
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
// -----------------------------------------------------------------------------
void HistogramPlotter::Draw_WeightLabels(TCanvas *c) {
    c->cd();

    TPaveText *pt, *pt2;

    pt = new TPaveText(0.840, 0.44, 0.91, 0.44, "NDC");
    pt->AddText(Form("#muB Genie Tune"));
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    
    if (_util.weight_tune)
        pt->Draw();

    pt2 = new TPaveText(0.839, 0.40, 0.909, 0.40, "NDC");
    pt2->AddText(Form("PPFX CV Corr."));
    pt2->SetBorderSize(0);
    pt2->SetFillColor(0);
    pt2->SetFillStyle(0);
    pt2->SetTextSize(0.03);
    
    if (_util.weight_ppfx)
        pt2->Draw();
}
// -----------------------------------------------------------------------------
void HistogramPlotter::Draw_VarMode(TCanvas *c ) {
    c->cd();

    TPaveText *pt;

    pt = new TPaveText(0.46, 0.88, 0.7, 0.95, "NDC");
    pt->AddText(_util.variation);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->SetTextColor(kViolet + 2);
    
    if (std::string(_util.variation) != "empty")
        pt->Draw();
}
// -----------------------------------------------------------------------------
void HistogramPlotter::Draw_Area_Norm(TCanvas *c) {
    c->cd();

    TPaveText *pt;

    pt = new TPaveText(0.3215, 0.936, 0.3215, 0.936, "NDC");
    pt->AddText("Area Normalised");
    pt->SetTextColor(kGreen + 2);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}
// -----------------------------------------------------------------------------
bool HistogramPlotter::GetHistograms(std::vector<TH1D *> &hist, std::string hist_name, std::string cut_name, std::string plotmode, bool &found_data, bool &found_ext, bool &found_dirt) {

    // Plots are by classification
    if (plotmode == "classifications") {

        // Loop over the classifications and get the histograms
        for (unsigned int i = 0; i < _util.classification_dirs.size(); i++) {

            // Data
            if (i == _util.k_leg_data) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
                
                if (hist.at(i) == NULL) {
                    found_data = false;
                }
            }
            // EXT
            else if (i == _util.k_leg_ext) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
                
                if (hist.at(i) == NULL){
                    found_ext = false;
                }
            }
            // Dirt
            else if (i == _util.k_leg_dirt) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));

                if (hist.at(i) == NULL) {
                    found_dirt = false;
                }
            }
            // MC
            else {

                // MC
                if (hist.at(i) != NULL && (i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt))
                    continue;

                _util.GetHist(f_nuexsec, hist.at(i), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));

                // Must have MC for this to work for now...
                if (hist.at(i) == NULL)
                    return false;
            }
        }
    }
    // By particles type
    else if (plotmode == "particle") {

        // Loop over the classifications and get the histograms
        for (unsigned int i = 0; i < _util.particle_types.size(); i++) {

            // Data
            if (i == _util.k_part_data) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("ParticleStack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.particle_types.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.particle_types.at(i).c_str()));
                
                if (hist.at(i) == NULL) {
                    found_data = false;
                }
            }
            // EXT
            else if (i == _util.k_part_ext) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("ParticleStack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.particle_types.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.particle_types.at(i).c_str()));
                
                if (hist.at(i) == NULL) {
                    found_ext = false;
                }
            }
            // Dirt
            else if (i == _util.k_part_dirt) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("ParticleStack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.particle_types.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.particle_types.at(i).c_str()));

                if (hist.at(i) == NULL) {
                    found_dirt = false;
                }
            }
            // MC
            else {

                // MC
                if (hist.at(i) != NULL && (i == _util.k_part_data || i == _util.k_part_ext || i == _util.k_part_dirt))
                    continue;

                _util.GetHist(f_nuexsec, hist.at(i), Form("ParticleStack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.particle_types.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.particle_types.at(i).c_str()));

                // Must have MC for this to work for now...
                if (hist.at(i) == NULL)
                    return false;
            }
        }
    }
    // Pi0 type
    else if (plotmode == "classifications_pi0") {

        // Loop over the classifications and get the histograms
        for (unsigned int i = 0; i < _util.classification_dirs.size(); i++) {
            // Data
            if (i == _util.k_leg_data){

                _util.GetHist(f_nuexsec, hist.at(i), Form("pizero/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));
                
                if (hist.at(i) == NULL){
                    found_data = false;
                }
            }
            // EXT
            else if (i == _util.k_leg_ext) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("pizero/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));
                
                if (hist.at(i) == NULL) {
                    found_ext = false;
                }
            }
            // Dirt
            else if (i == _util.k_leg_dirt) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("pizero/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));

                if (hist.at(i) == NULL) {
                    found_dirt = false;
                }
            }
            // MC
            else {

                // MC
                if (hist.at(i) != NULL && (i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt))
                    continue;

                _util.GetHist(f_nuexsec, hist.at(i), Form("pizero/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));

                // Must have MC for this to work for now...
                if (hist.at(i) == NULL)
                    return false;
            }
        
        }

    }
    // Pi0 type
    else {

        // Loop over the classifications and get the histograms
        for (unsigned int i = 0; i < _util.classification_dirs.size(); i++) {
            // Data
            if (i == _util.k_leg_data) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("numu/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));
                if (hist.at(i) == NULL) {
                    found_data = false;
                }
            }
            // EXT
            else if (i == _util.k_leg_ext){

                _util.GetHist(f_nuexsec, hist.at(i), Form("numu/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));
                
                if (hist.at(i) == NULL) {
                    found_ext = false;
                }
            }
            // Dirt
            else if (i == _util.k_leg_dirt) {

                _util.GetHist(f_nuexsec, hist.at(i), Form("numu/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));

                if (hist.at(i) == NULL) {
                    found_dirt = false;
                }
            }
            // MC
            else {

                // MC
                if (hist.at(i) != NULL && (i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt))
                    continue;

                _util.GetHist(f_nuexsec, hist.at(i), Form("numu/%s_%s", hist_name.c_str(), _util.classification_dirs.at(i).c_str()));

                // Must have MC for this to work for now...
                if (hist.at(i) == NULL)
                    return false;
            }
        
        }

    }

    return true;
}
// -----------------------------------------------------------------------------
void HistogramPlotter::SetFillColours(std::vector<TH1D *> &hist, std::string plotmode, bool found_data, bool found_dirt, bool found_ext, unsigned int k_plot_data, unsigned int k_plot_ext, unsigned int k_plot_dirt) {

    if (plotmode == "classifications" || plotmode == "classifications_pi0"|| plotmode == "classifications_numu") {
        hist.at(_util.k_nue_cc)->SetFillColor(30);
        hist.at(_util.k_nuebar_cc)->SetFillColor(32);
        hist.at(_util.k_numu_cc)->SetFillColor(28);
        hist.at(_util.k_nc_pi0)->SetFillColor(36);
        hist.at(_util.k_cosmic)->SetFillColor(1);
        hist.at(_util.k_nc)->SetFillColor(46);
        hist.at(_util.k_nu_out_fv)->SetFillColor(kViolet - 7);
        hist.at(_util.k_numu_cc_pi0)->SetFillColor(42);
        hist.at(_util.k_unmatched)->SetFillColor(12);
    }
    // By particle type
    else {

        hist.at(_util.k_electron)->SetFillColor(30);
        hist.at(_util.k_neutron)->SetFillColor(38);
        hist.at(_util.k_pion)->SetFillColor(28);
        hist.at(_util.k_photon)->SetFillColor(42);
        hist.at(_util.k_muon)->SetFillColor(36);
        hist.at(_util.k_part_cosmic)->SetFillColor(1);
        hist.at(_util.k_proton)->SetFillColor(46);
        hist.at(_util.k_kaon)->SetFillColor(kViolet - 7);
        hist.at(_util.k_part_unmatched)->SetFillColor(12);
    }

    // This is generic
    if (found_data) {
        hist.at(k_plot_data)->SetMarkerStyle(20);
        hist.at(k_plot_data)->SetMarkerSize(0.5);
    }
    if (found_ext) {
        hist.at(k_plot_ext)->SetFillColor(41);
        hist.at(k_plot_ext)->SetFillStyle(3345);
    }

    if (found_dirt) {
        hist.at(k_plot_dirt)->SetFillColor(2);
        hist.at(k_plot_dirt)->SetFillStyle(3354);
    }
}
// -----------------------------------------------------------------------------
void HistogramPlotter::SetLegend(std::vector<TH1D *> hist, TLegend *leg_stack, std::vector<double> hist_integrals, bool found_data, bool found_dirt, bool found_ext, unsigned int k_plot_data, unsigned int k_plot_ext, unsigned int k_plot_dirt, std::string plotmode){

    if (plotmode == "classifications" || plotmode == "classifications_pi0" || plotmode == "classifications_numu") {
        if (found_data) leg_stack->AddEntry(hist.at(k_plot_data), Form("Beam-On Data (%2.1f)", hist_integrals.at(_util.k_leg_data)), "lep");
        if (found_dirt) leg_stack->AddEntry(hist.at(k_plot_dirt), Form("Out-of Cryo (%2.1f)", hist_integrals.at(_util.k_leg_dirt)), "f");
        if (found_ext)  leg_stack->AddEntry(hist.at(k_plot_ext),  Form("Beam-Off Data (%2.1f)", hist_integrals.at(_util.k_leg_ext)), "f");
        // leg_stack->AddEntry(hist.at(_util.k_unmatched),       Form("Unmatched (%2.1f)",           hist_integrals.at(_util.k_unmatched)),    "f"); // This should be zero, so dont plot
        leg_stack->AddEntry(hist.at(_util.k_nc_pi0),      Form("NC #pi^{0} (%2.1f)", hist_integrals.at(_util.k_nc_pi0)), "f");
        leg_stack->AddEntry(hist.at(_util.k_nc),          Form("NC (%2.1f)", hist_integrals.at(_util.k_nc)), "f");
        leg_stack->AddEntry(hist.at(_util.k_numu_cc_pi0), Form("#nu_{#mu} CC #pi^{0} (%2.1f)", hist_integrals.at(_util.k_numu_cc_pi0)), "f");
        leg_stack->AddEntry(hist.at(_util.k_numu_cc),     Form("#nu_{#mu} CC (%2.1f)", hist_integrals.at(_util.k_numu_cc)), "f");
        leg_stack->AddEntry(hist.at(_util.k_cosmic),      Form("Cosmic (%2.1f)", hist_integrals.at(_util.k_cosmic)), "f");
        leg_stack->AddEntry(hist.at(_util.k_nu_out_fv),   Form("#nu OutFV (%2.1f)", hist_integrals.at(_util.k_nu_out_fv)), "f");
        leg_stack->AddEntry(hist.at(_util.k_nuebar_cc),   Form("#bar{#nu}_{e} CC (%2.1f)", hist_integrals.at(_util.k_nuebar_cc)), "f");
        leg_stack->AddEntry(hist.at(_util.k_nue_cc),      Form("#nu_{e} CC (%2.1f)", hist_integrals.at(_util.k_nue_cc)), "f");
    }

    else {
        if (found_data) leg_stack->AddEntry(hist.at(k_plot_data), Form("Beam-On Data (%2.1f)", hist_integrals.at(_util.k_part_data)), "lep");
        if (found_dirt) leg_stack->AddEntry(hist.at(k_plot_dirt), Form("Out-of Cryo (%2.1f)", hist_integrals.at(_util.k_part_dirt)), "f");
        if (found_ext)  leg_stack->AddEntry(hist.at(k_plot_ext),  Form("Beam-Off Data (%2.1f)", hist_integrals.at(_util.k_part_ext)), "f");
        // leg_stack->AddEntry(hist.at(_util.k_part_unmatched), Form("Unmatched (%2.1f)", hist_integrals.at(_util.k_part_unmatched)), "f");
        leg_stack->AddEntry(hist.at(_util.k_kaon),        Form("K (%2.1f)", hist_integrals.at(_util.k_kaon)), "f");
        leg_stack->AddEntry(hist.at(_util.k_proton),      Form("p (%2.1f)", hist_integrals.at(_util.k_proton)), "f");
        leg_stack->AddEntry(hist.at(_util.k_part_cosmic), Form("Cosmic (%2.1f)", hist_integrals.at(_util.k_part_cosmic)), "f");
        leg_stack->AddEntry(hist.at(_util.k_muon),        Form("#mu (%2.1f)", hist_integrals.at(_util.k_muon)), "f");
        leg_stack->AddEntry(hist.at(_util.k_photon),      Form("#gamma (%2.1f)", hist_integrals.at(_util.k_photon)), "f");
        leg_stack->AddEntry(hist.at(_util.k_pion),        Form("#pi (%2.1f)", hist_integrals.at(_util.k_pion)), "f");
        leg_stack->AddEntry(hist.at(_util.k_neutron),     Form("n (%2.1f)", hist_integrals.at(_util.k_neutron)), "f");
        leg_stack->AddEntry(hist.at(_util.k_electron),    Form("e (%2.1f)", hist_integrals.at(_util.k_electron)), "f");
    }
}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeStack(std::string hist_name, std::string cut_name, bool _area_norm, bool logy, double y_scale_factor, const char *x_axis_name,
                                  const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, double Data_POT,
                                  const char *print_name, bool override_data_mc_comparison, std::string plotmode, bool plotvar,
                                  bool centerxaxis, bool scale) {

    // Check if we want this plot for variations
    if (!plotvar && std::string(_util.variation) != "empty") return;

    std::vector<TH1D *> hist;           // The vector of histograms from the file for the plot
    std::vector<double> hist_integrals; // The integrals of all the histograms

    unsigned int k_plot_data, k_plot_ext, k_plot_dirt; // Maps to the correct enum value between the classifications and particle type

    if (plotmode == "classifications" || plotmode == "classifications_pi0" || plotmode == "classifications_numu") {
        hist.resize(_util.k_classifications_MAX);
        hist_integrals.resize(_util.k_classifications_MAX, 0.0);
        k_plot_data = _util.k_leg_data;
        k_plot_ext = _util.k_leg_ext;
        k_plot_dirt = _util.k_leg_dirt;
    }
    else {
        hist.resize(_util.k_particles_MAX);
        hist_integrals.resize(_util.k_particles_MAX, 0.0);
        k_plot_data = _util.k_part_data;
        k_plot_ext = _util.k_part_ext;
        k_plot_dirt = _util.k_part_dirt;
    }

    double integral_mc_ext{0.0}; // Integral of MC + EXT -- needs to be removed
    double y_maximum{0};         // y max for scaling histogram scale

    std::vector<double> chi2;
    TPaveText *pt_bottom;

    TH1D *h_ratio;
    TH1D *h_ratio_error;
    TH1D *h_mc_ext_sum;

    // bools for checking if plots exist in the file
    bool found_data = true;
    bool found_ext = true;
    bool found_dirt = true;

    const bool p_value = false; // Choose whether to use a pvalue

    TPad *topPad;
    TPad *bottomPad;
    // TCanvas *c  = new TCanvas();
    TCanvas * c = new TCanvas(Form("c_%s_%s", print_name, cut_name.c_str()), "c", 500, 500);
    THStack *h_stack = new THStack();

    // Get the histograms from the file
    bool bool_hists = GetHistograms(hist, hist_name, cut_name, plotmode, found_data, found_ext, found_dirt);
    if (!bool_hists) {
        std::cout << "failed to get the histograms! " << hist_name << std::endl;
        return;
    }

    // If true, will turn off plotting of data and off beam
    if (override_data_mc_comparison) {
        found_data = false;
        found_ext = false;
    }

    // If there is data and ext, then we make the ratio plot too
    if (found_data && found_ext) {
        topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

        _util.SetTPadOptions(topPad, bottomPad);
    }
    // Otherwise just use an unsplit canvas
    else {
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.20);
    }

    // Scaling and setting stats
    for (unsigned int i = 0; i < hist.size(); i++) {

        if (i == k_plot_data) {
            if (found_data) {
                hist.at(i)->SetStats(kFALSE);
                hist_integrals.at(i) = hist.at(i)->Integral();
            }
        }

        // Scale EXT
        else if (i == k_plot_ext) {
            if (found_ext) {

                hist.at(i)->SetStats(kFALSE);
                if (scale) hist.at(i)->Scale(_util.ext_scale_factor);
                hist_integrals.at(i) = hist.at(i)->Integral();
                integral_mc_ext += hist.at(i)->Integral();
            }
        }

        // Scale Dirt
        else if (i == k_plot_dirt){
            if (found_dirt) {

                hist.at(i)->SetStats(kFALSE);
                if (scale) hist.at(i)->Scale(_util.dirt_scale_factor);
                hist_integrals.at(i) = hist.at(i)->Integral();
                integral_mc_ext += hist.at(i)->Integral();
            }
        }

        // Scale MC
        else {

            hist.at(i)->SetStats(kFALSE);
            if (scale) hist.at(i)->Scale(_util.mc_scale_factor);
            hist_integrals.at(i) = hist.at(i)->Integral();
            integral_mc_ext += hist.at(i)->Integral();
            // std::cout <<  hist.at(i)->Integral()<< std::endl;
        }
    }

    // Set fill colours of stacked histogram
    SetFillColours(hist, plotmode, found_data, found_dirt, found_ext, k_plot_data, k_plot_ext, k_plot_dirt);

    // Normalisation by area
    if (_area_norm && found_data) {

        if (integral_mc_ext != 0) {

            for (unsigned int i = 0; i < hist.size(); i++) {
                if (i == k_plot_data)
                    continue; // Dont scale the data
                if (i == 0)
                    std::cout << "area norm scale factor: " << hist_integrals.at(k_plot_data) / integral_mc_ext << std::endl;
                hist.at(i)->Scale(hist_integrals.at(k_plot_data) / integral_mc_ext);
            }
        }
    }

    // Add the histograms to the stack
    for (unsigned int i = 0; i < hist.size(); i++) {
        if (i == k_plot_data)
            continue; // Dont use the data
        if (i == k_plot_ext && !found_ext)
            continue; // Skip ext if not there
        if (i == k_plot_dirt && !found_dirt)
            continue; // skip dirt if not there

        h_stack->Add(hist.at(i));
    }

    // Get the maximum y value for scaling
    if (found_data)
        y_maximum = std::max(hist.at(k_plot_data)->GetMaximum(), h_stack->GetMaximum());
    else
        y_maximum = h_stack->GetMaximum();

    // Set the axis in the case of a log plot
    if (logy == true && found_data) {

        if (hist.at(0)->GetMinimum() != 0.0)
            hist.at(k_plot_data)->SetMinimum(hist.at(0)->GetMinimum() / 2.); // hist.at(0) is just checking the bounds of the first histogram

        h_stack->SetMinimum(2.5);
        
        if (hist.at(0)->GetMinimum() == 0.0)
            hist.at(k_plot_data)->SetMinimum(hist.at(0)->GetMinimum() + 0.0001 / 2.);

        h_stack->SetMaximum(y_maximum * y_scale_factor);

        h_stack->Draw("hist");
    }
    else if (logy && !found_data) {
        h_stack->SetMinimum(1);
        c->SetLogy();
        h_stack->Draw("hist");
    }
    else { // Set the axis in the case of a non-log plot
        h_stack->SetMinimum(0);
        h_stack->SetMaximum(y_maximum * y_scale_factor);
        h_stack->Draw("hist");
    }

    // Set the y axis of the stack
    if (!_area_norm)
        h_stack->GetYaxis()->SetTitle("Entries");
    else
        h_stack->GetYaxis()->SetTitle("Entries [A.U.]");

    // Customise the stacked histogram
    if (!found_data){
        h_stack->GetYaxis()->SetTitleSize(0.04);
        h_stack->GetYaxis()->SetLabelSize(0.04);
    }
    else {
        h_stack->GetYaxis()->SetTitleSize(0.05);
        h_stack->GetYaxis()->SetLabelSize(0.05);
    }
    
    if (found_data)
        h_stack->GetXaxis()->SetLabelSize(0);

    // Customise bin labels for pass and fail type of histograms (ones with 2 bins)
    if (h_stack->GetXaxis()->GetNbins() == 2) {
        h_stack->GetXaxis()->SetNdivisions(2, 0, 0, kFALSE);
        h_stack->GetXaxis()->CenterLabels(kTRUE);
        h_stack->GetXaxis()->SetBinLabel(1, "Fail");
        h_stack->GetXaxis()->SetBinLabel(2, "Pass");
    }

    // Set Poisson Errors for data
    // hist.at(_util.k_leg_data)->SetBinErrorOption(TH1::kPoisson);

    if (found_data && found_ext)
        h_stack->GetXaxis()->SetLabelOffset(10);
    else
        h_stack->GetXaxis()->SetTitle(x_axis_name);

    if (found_data)
        hist.at(k_plot_data)->Draw("same PE");

    // MC error histogram ------------------------------------------------------
    TH1D *h_error_hist = (TH1D *)hist.at(0)->Clone("h_error_hist");

    for (unsigned int i = 0; i < hist.size(); i++) {
        if (i == k_plot_data || i == 0)
            continue; // Dont use the data
        if (i == k_plot_ext && !found_ext)
            continue; // Skip ext if not there
        if (i == k_plot_dirt && !found_dirt)
            continue; // skip dirt if not there

        if (i == 0)
            continue; // Aleady got this histogram from the clone

        h_error_hist->Add(hist.at(i), 1);
    }

    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->SetLineWidth(0);
    h_error_hist->Draw("e2, same");
    if (hist_name == "h_reco_track_multiplicity")h_error_hist->Draw("e2, same, TEXT00"); // for the track multiplicity plot draw the event count so we can get the % of 0 track events


    // Set the legend ----------------------------------------------------------
    TLegend *leg_stack;
    if (found_data){
        leg_stack = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
        leg_stack->SetNColumns(2);
    }
    else {
        leg_stack = new TLegend(0.8, 0.87, 0.98, 0.32);
    }
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    SetLegend(hist, leg_stack, hist_integrals, found_data, found_dirt, found_ext, k_plot_data, k_plot_ext, k_plot_dirt, plotmode);

    leg_stack->Draw();

    if (!logy)
        h_stack->GetYaxis()->SetRangeUser(0, y_maximum * y_scale_factor);
    else if (logy && found_data)
        topPad->SetLogy();

    // Calculate the chi2
    if (found_data) {
        TH1D *h_last = (TH1D *)h_stack->GetStack()->Last();
        chi2 = Chi2Calc(h_last, hist.at(k_plot_data), _area_norm, hist_integrals.at(k_plot_data));
    }

    // Now create the ratio of data to MC ----------------------------------
    if (found_data) {

        bottomPad->cd();

        h_ratio = (TH1D *)hist.at(k_plot_data)->Clone("h_ratio");
        h_mc_ext_sum = (TH1D *)hist.at(0)->Clone("h_mc_ext_sum");

        for (unsigned int i = 0; i < hist.size(); i++) {
            if (i == k_plot_data || i == 0)
                continue; // Dont use the data and nue cc because already been cloned
            h_mc_ext_sum->Add(hist.at(i), 1);
        }

        // h_ratio->Add(h_mc_ext_sum, -1); // Turn off for data / MC + ext
        h_ratio->Divide(h_mc_ext_sum);

        h_ratio->GetXaxis()->SetLabelSize(0.13);
        h_ratio->GetYaxis()->SetLabelSize(0.13);
        h_ratio->GetXaxis()->SetTitleOffset(0.9);
        h_ratio->GetXaxis()->SetTitleSize(0.13);

        // Customise bin labels for pass and fail type of histograms (ones with 2 bins)
        if (h_ratio->GetXaxis()->GetNbins() == 2) {
            h_ratio->GetXaxis()->SetNdivisions(2, 0, 0, kFALSE);
            h_ratio->GetXaxis()->CenterLabels(kTRUE);
            h_ratio->GetXaxis()->SetBinLabel(1, "Fail");
            h_ratio->GetXaxis()->SetBinLabel(2, "Pass");

            // Its a veto in this case so swap them
            if (hist_name == "h_reco_crtveto") {
                h_ratio->GetXaxis()->SetBinLabel(1, "Pass");
                h_ratio->GetXaxis()->SetBinLabel(2, "Fail");
            }
        }

        h_ratio->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

        // For percent difference
        // h_ratio->GetYaxis()->SetTitle("(Data - MC) / MC ");
        // h_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);
        // h_ratio->GetYaxis()->SetRangeUser(-0.2, 0.2);

        // For ratio
        if (cut_name == "Unselected" || cut_name == "SoftwareTrig" || cut_name == "Slice_ID" || cut_name == "Topo_Score" || cut_name == "In_FV" || cut_name == "e_candidate")
            h_ratio->GetYaxis()->SetRangeUser(0.8, 1.2);
        else
            h_ratio->GetYaxis()->SetRangeUser(0, 2.0);

        h_ratio->GetYaxis()->SetTitle("#frac{Beam-On}{(Overlay + Beam-Off)}");

        h_ratio->GetXaxis()->SetTitle(x_axis_name);
        h_ratio->GetYaxis()->SetTitleSize(10);
        h_ratio->GetYaxis()->SetTitleFont(44);
        h_ratio->GetYaxis()->CenterTitle();
        h_ratio->GetYaxis()->SetTitleOffset(2.5);
        h_ratio->SetTitle(" ");
        // h_ratio->SetBinErrorOption(TH1::kPoisson);

        h_ratio->Draw("E0,same");

        // Draw the error hist
        h_ratio_error = (TH1D *)h_error_hist->Clone("h_ratio_error");
        h_ratio_error->Divide(h_ratio_error);
        h_ratio_error->Draw("e2, same");
    
        // Choose whether to center the xaxis labels. Only makes sense for counting type of plots
        if (centerxaxis) {
            h_ratio->GetXaxis()->CenterLabels(kTRUE);
        }

        // Now doing this stuff on the bottom pad
        //x_min, y_min, x_max, y_max
        // Reduced chi2
        pt_bottom = new TPaveText(.12, .80, .30, .96, "NBNDC");
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
    }

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Draw variation mode if in that mode
    Draw_VarMode(c);

    // Draw other data specifc quantities
    if (found_data) {

        // Turn this off for now until we understand this
        // _util.Draw_Data_MC_Ratio(c, hist_integrals.at(_util.k_leg_data)/integral_mc_ext, 0.34, 0.936, 0.34, 0.936 );

        _util.Draw_Data_POT(c, Data_POT, 0.45, 0.915, 0.45, 0.915);

        if (_area_norm)
            Draw_Area_Norm(c);

        // Add the weight labels
        // Draw_WeightLabels(c);
    }

    if (plotmode == "classifications_pi0" || plotmode == "classifications_numu"){
        c->Print(print_name);
        return;
    }

    if (std::string(_util.variation) == "empty")
        c->Print(Form("plots/run%s/%s", _util.run_period, print_name));
    else
        c->Print(Form("plots/run%s/detvar/%s/%s", _util.run_period, _util.variation, print_name));

    // delete c;
}
// -----------------------------------------------------------------------------
void HistogramPlotter::CallMakeStack(int cut_index, double Data_POT) {

    // MakeStack(std::string hist_name, std::string cut_name, bool area_norm, bool logy, const char* x_axis_name, double y_scale_factor,
    //                             const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, const char* print_name, bool override_data_mc_comparison )

    // Reco X SCE
    MakeStack("h_reco_vtx_x_sce", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.3, "Reco Vertex X (Space Charge Corr) [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_vtx_x_sce.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Reco Y SCE
    MakeStack("h_reco_vtx_y_sce", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Reco Vertex Y (Space Charge Corr) [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_vtx_y_sce.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Reco Z SCE
    MakeStack("h_reco_vtx_z_sce", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.2, "Reco Vertex Z (Space Charge Corr) [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_vtx_z_sce.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Leading Shower Momentum
    MakeStack("h_reco_leading_mom", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower Momentum [GeV/c]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_leading_mom.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // 2D distance shower vertex to reco nu vertex
    MakeStack("h_reco_shower_to_vtx_dist", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Shower to Vertex Distance [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_to_vtx_dist.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Leading Shower hits in all planes
    MakeStack("h_reco_shr_hits_max", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower Hits All Planes",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_hits_max.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Leading shower phi
    MakeStack("h_reco_leading_shower_phi", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower Phi [degrees]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_leading_shower_phi.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Leading shower theta
    MakeStack("h_reco_leading_shower_theta", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower Theta [degrees]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_leading_shower_theta.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Leading shower cos theta
    MakeStack("h_reco_leading_shower_cos_theta", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.5, "Leading Shower Cos(#theta)",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_leading_shower_cos_theta.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Leading shower multiplicity
    MakeStack("h_reco_shower_multiplicity", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.5, "Shower Multiplicty",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_multiplicity.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, true, true);
    MakeStack("h_reco_shower_multiplicity", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 80, "Shower Multiplicty",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_multiplicity_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, true, false);

    // Leading track multiplicity
    MakeStack("h_reco_track_multiplicity", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.5, "Track Multiplicty",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_track_multiplicity.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, true, true);

    // Topological Score
    MakeStack("h_reco_topological_score", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Topological Score",  0.3, 0.8, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_topological_score.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);
    MakeStack("h_reco_topological_score", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 50, "Topological Score",  0.3, 0.8, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_topological_score_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, false);

    // Track shower dist
    MakeStack("h_reco_track_shower_dist", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Longest Track Leading Shower Distance [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_track_shower_dist.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Track shower angle
    MakeStack("h_reco_track_shower_angle", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.5, "Longest Track Leading Shower Angle [degrees]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_track_shower_angle.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Ratio hits from showers to slice
    MakeStack("h_reco_hits_ratio", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Hit Ratio",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_hits_ratio.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Shower score
    MakeStack("h_reco_shower_score", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Shower Score",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_score.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    MakeStack("h_reco_shower_score", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 1.0, "Shower Score",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_score_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, false);

    // Track score
    MakeStack("h_reco_track_score", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Track Score",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_track_score.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Calibrated energy of all the showers
    MakeStack("h_reco_shower_energy_tot_cali", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Total Calibrated Energy of all Showers [GeV]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_energy_tot_cali.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    
    // Calibrated energy of all the showers with uneaven bins
    MakeStack("h_reco_shower_energy_tot_cali_rebin", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Total Calibrated Energy of all Showers [GeV]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_energy_tot_cali_rebin.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Calibrated energy of the leading shower
    MakeStack("h_reco_shower_energy_cali", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Reco Energy Leading Shower [GeV]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_energy_cali.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    
    // Calibrated energy of all the showers with uneaven bins
    MakeStack("h_reco_shower_energy_cali_rebin", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.6, "Reco Energy Leading Shower [GeV]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shower_energy_cali_rebin.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);


    // Total number of hits for all showers
    MakeStack("h_reco_shr_hits_tot", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Total Num of hits for all Showers",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_hits_tot.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Total number of hits for the all showers in the collection plane
    MakeStack("h_reco_shr_hits_y_tot", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Total Num of hits for all Showers in Collection Plane",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_hits_y_tot.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Optical Filter Beam
    MakeStack("h_reco_opfilter_beam", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Common Optical Filter PE [PE]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_opfilter_beam.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    MakeStack("h_reco_opfilter_beam", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 1.0, "Common Optical Filter PE [PE]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_opfilter_beam_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, false);

    // Optical Filter Veto
    MakeStack("h_reco_opfilter_veto", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Common Optical Filter Michel Veto [PE]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_opfilter_veto.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    MakeStack("h_reco_opfilter_veto", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 1.0, "Common Optical Filter Michel Veto [PE]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_opfilter_veto_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, false);

    // Software trigger
    MakeStack("h_reco_softwaretrig", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Software Trigger",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_softwaretrig.pdf", _util.cut_dirs.at(cut_index).c_str()), true, "classifications", false, true, true);
    // Software trigger
    MakeStack("h_reco_softwaretrig", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 1.0, "Software Trigger",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_softwaretrig_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), true, "classifications", false, true, false);

    // Slice ID
    MakeStack("h_reco_nslice", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.6, "Pandora Slice ID",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_nslice.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    MakeStack("h_reco_nslice", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 800, "Pandora Slice ID",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_nslice_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, false);

    // Slice Cluster Fraction
    MakeStack("h_reco_slclustfrac", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Slice Cluster Fraction",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_slclustfrac.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Cosmic Inpact Parameter in 3D
    MakeStack("h_reco_CosmicIPAll3D", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.7, "Pandora Cosmic Impact Parameter 3D [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_CosmicIPAll3D.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // 3D direction between closest pfp not in the neutrino slice
    MakeStack("h_reco_CosmicDirAll3D", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "3D Pandora Cosmic Impact Direction",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_CosmicDirAll3D.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Shower dEdx with the track fitter u plane
    MakeStack("h_reco_shr_tkfit_dedx_u", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "U Plane dEdx (track fitter) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_u.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Shower dEdx with the track fitter v plane
    MakeStack("h_reco_shr_tkfit_dedx_v", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "V Plane dEdx (track fitter) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_v.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Shower dEdx with the track fitter y plane
    MakeStack("h_reco_shr_tkfit_dedx_y", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower dEdx (Collection Plane) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_y.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Shower dEdx with the track fitter plane with the most hits
    MakeStack("h_reco_shr_tkfit_dedx_max", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower dEdx (All Planes) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_max.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Shower dEdx with the track fitter plane with the most hits for events with tracks only
    MakeStack("h_reco_shr_tkfit_dedx_max_with_tracks", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower dEdx (All Planes) (with tracks) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_max_with_tracks.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Shower dEdx with the track fitter y plane when there is no tracks
    MakeStack("h_reco_shr_tkfit_dedx_y_no_tracks", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower dEdx (Collection Plane) (0 tracks) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_y_no_tracks.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Shower dEdx with the track fitter plane with most hits when there is no tracks
    MakeStack("h_reco_shr_tkfit_dedx_max_no_tracks", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.5, "Leading Shower dEdx (All Planes) (0 tracks) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_max_no_tracks.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Shower Flash Time
    MakeStack("h_reco_flash_time", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Flash Time [#mus]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_flash_time.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Shower Flash PE
    MakeStack("h_reco_flash_pe", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Largest Flash Intensity [PE]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_flash_pe.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Shower Subcluster All Planes
    MakeStack("h_reco_shrsubclusters", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower Sub Cluster All Planes",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shrsubclusters.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Shower Moliere Average
    MakeStack("h_reco_shrmoliereavg", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.6, "Leading Shower Moliere Average [deg]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shrmoliereavg.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Median of 1st component of shr PCA (5cm window)
    MakeStack("h_reco_shrMCSMom", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower MCS Momentum [MeV/c]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shrMCSMom.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Distance from closest CRT tagged cosmic to neutrino vertex
    MakeStack("h_reco_closestNuCosmicDist", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Closest CRT Tagged Cosmic to #nu Vertex Distance [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_closestNuCosmicDist.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Longest Track Length
    MakeStack("h_reco_trk_len", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Longest Track Length [cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_trk_len.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Reco Neutrino Energy
    MakeStack("h_reco_nu_e", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Reconstructed Neutrino Energy [GeV]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_nu_e.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", true, false, true);

    // Contained Fraction
    MakeStack("h_reco_contained_fraction", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Contained Fraction (PFP hits in FV / hits in slice)",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_contained_fraction.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    MakeStack("h_reco_contained_fraction", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 1.0, "Contained Fraction (PFP hits in FV / hits in slice)",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_contained_fraction_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, false);

    // Run number
    MakeStack("h_reco_run_number", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Run Number",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_run_number.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // Neutrino Purity
    MakeStack("h_reco_nu_purity_from_pfp", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Neutrino Purity",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_nu_purity_from_pfp.pdf", _util.cut_dirs.at(cut_index).c_str()), true, "classifications", false, false, true);
    MakeStack("h_reco_nu_purity_from_pfp", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 1.0, "Neutrino Purity",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_nu_purity_from_pfp_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), true, "classifications", false, false, false);

    // CRT Veto
    MakeStack("h_reco_crtveto", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "CRT Veto",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_crtveto.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);

    // CRT HitPE
    MakeStack("h_reco_crthitpe", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "CRT Hit Intensity [PE]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_crthitpe.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, true);
    MakeStack("h_reco_crthitpe", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, true, 1.0, "CRT Hit Intensity [PE]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_crthitpe_logy.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "classifications", false, false, false);

    // Stacked Histograms by particle type

    // dEdx Y Plane using trackfits
    MakeStack("h_reco_shr_tkfit_dedx_y_par", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower dEdx (Collection Plane) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_y_par.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "particle", true, false, true);
    
    // dEdx max Plane using trackfits
    MakeStack("h_reco_shr_tkfit_dedx_max_par", _util.cut_dirs.at(cut_index).c_str(),
              area_norm, false, 1.0, "Leading Shower dEdx (All Planes) [MeV/cm]",  0.35, 0.85, 0.55, 0.85, Data_POT,
              Form("cuts/%s/reco_shr_tkfit_dedx_max_par.pdf", _util.cut_dirs.at(cut_index).c_str()), false, "particle", true, false, true);


}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeFlashPlot(double Data_POT, const char *print_name, std::string histname) {

    std::vector<TH1D *> hist(_util.k_type_MAX);
    std::vector<double> hist_integrals(_util.k_type_MAX, 0.0); // The integrals of all the histograms
    double integral_mc_ext = 0.0;

    TH1D *h_ratio;
    TH1D *h_ratio_error;
    TH1D *h_mc_ext_sum;

    TPad *topPad;
    TPad *bottomPad;
    TCanvas * c = new TCanvas(Form("c_%s", histname.c_str()), "c", 500, 500);
    THStack *h_stack = new THStack();

    for (unsigned int k = 0; k < _util.type_prefix.size(); k++) {
        _util.GetHist(f_nuexsec, hist.at(k), Form("Flash/%s_%s", histname.c_str(), _util.type_prefix.at(k).c_str()));
        
        if (hist.at(k) == NULL) {
            std::cout << "Couldn't get all the flash histograms so exiting function..." << std::endl;
            return;
        }
    }

    topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    _util.SetTPadOptions(topPad, bottomPad);

    for (unsigned int i = 0; i < hist.size(); i++) {

        if (i == _util.k_data) {

            hist.at(i)->SetStats(kFALSE);
            hist.at(_util.k_data)->SetMarkerStyle(20);
            hist.at(_util.k_data)->SetMarkerSize(0.5);
            hist_integrals.at(_util.k_data) = hist.at(_util.k_data)->Integral();
        }

        // Scale EXT
        else if (i == _util.k_ext) {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(_util.ext_scale_factor);
            hist.at(_util.k_ext)->SetFillColor(41);
            hist.at(_util.k_ext)->SetFillStyle(3345);
            hist_integrals.at(_util.k_ext) = hist.at(_util.k_ext)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_ext);
        }

        // Scale Dirt
        else if (i == _util.k_dirt) {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(_util.dirt_scale_factor);
            hist.at(_util.k_dirt)->SetFillColor(2);
            hist.at(_util.k_dirt)->SetFillStyle(3354);
            hist_integrals.at(_util.k_dirt) = hist.at(_util.k_dirt)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_dirt);
        }

        // Scale MC
        else {
            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->Scale(_util.mc_scale_factor);
            hist.at(i)->SetFillColor(30);
            hist_integrals.at(_util.k_mc) = hist.at(_util.k_mc)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_mc);
        }
    }

    // Normalisation by area
    if (area_norm) {

        if (integral_mc_ext != 0) {

            for (unsigned int i = 0; i < hist.size(); i++) {
                if (i == _util.k_data)
                    continue; // Dont scale the data
                // if (i == 0) std::cout << "area norm scale factor: "  << hist_integrals.at(k_plot_data) / integral_mc_ext << std::endl;
                hist.at(i)->Scale(hist_integrals.at(_util.k_data) / integral_mc_ext);
            }
        }
    }

    // Add the histograms to the stack
    h_stack->Add(hist.at(_util.k_ext));
    h_stack->Add(hist.at(_util.k_mc));
    h_stack->Add(hist.at(_util.k_dirt));

    h_stack->Draw("hist");
    hist.at(_util.k_data)->Draw("same PE");

    if (area_norm)
        h_stack->GetYaxis()->SetTitle("Entries A.U. ");
    else
        h_stack->GetYaxis()->SetTitle("Entries");

    h_stack->GetYaxis()->SetTitleSize(0.05);
    h_stack->GetYaxis()->SetLabelSize(0.05);
    h_stack->GetXaxis()->SetLabelSize(0);
    h_stack->GetXaxis()->SetRangeUser(0, 23);
    if (histname == "h_flash_time_single_bin") h_stack->GetXaxis()->SetRangeUser(5.6,15.4);

    // MC error histogram ------------------------------------------------------
    TH1D *h_error_hist = (TH1D *)hist.at(_util.k_mc)->Clone("h_error_hist");

    for (unsigned int i = 0; i < hist.size(); i++) {
        if (i == _util.k_data)
            continue; // Dont use the data
        if (i == _util.k_mc)
            continue; // Aleady got this histogram from the clone

        h_error_hist->Add(hist.at(i), 1);
    }

    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2, same");

    TLegend *leg_stack;
    if (histname == "h_flash_time_single_bin")
        leg_stack = new TLegend(0.7, 0.2, 0.9, 0.45);
    else
        leg_stack = new TLegend(0.7, 0.6, 0.9, 0.85);
    
    leg_stack->SetBorderSize(0);
    
    if (histname != "h_flash_time_single_bin")
        leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(hist.at(_util.k_data), "Beam-On Data", "lep");
    leg_stack->AddEntry(hist.at(_util.k_dirt), "Out-of Cryo", "f");
    leg_stack->AddEntry(hist.at(_util.k_mc),   "Overlay", "f");
    leg_stack->AddEntry(hist.at(_util.k_ext),  "Beam-Off Data", "f");

    leg_stack->Draw();

    bottomPad->cd();

    h_ratio = (TH1D *)hist.at(_util.k_data)->Clone("h_ratio");
    h_mc_ext_sum = (TH1D *)hist.at(_util.k_mc)->Clone("h_mc_ext_sum");

    for (unsigned int i = 0; i < hist.size(); i++) {
        
        if (i == _util.k_data || i == _util.k_mc)
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
    h_ratio->GetXaxis()->SetRangeUser(0, 23);
    if (histname == "h_flash_time_single_bin") h_ratio->GetXaxis()->SetRangeUser(5.6,15.4);
    h_ratio->GetYaxis()->SetTitle("#frac{Beam-On}{(Overlay + Beam-Off)}");

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
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Draw Data to MC ratio
    _util.Draw_Data_MC_Ratio(c, double(hist_integrals.at(_util.k_data) * 1.0 / integral_mc_ext * 1.0), 0.34, 0.936, 0.34, 0.936);

    // Draw other data specifc quantities
    _util.Draw_Data_POT(c, Data_POT, 0.45, 0.915, 0.45, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);

    if (area_norm)
        Draw_Area_Norm(c);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeFlashPlotOMO(double Data_POT, const char *print_name, std::string histname) {

    std::vector<TH1D *> hist(_util.k_type_MAX);
    std::vector<double> hist_integrals(_util.k_type_MAX, 0.0); // The integrals of all the histograms
    double integral_mc_ext = 0.0;

    TH1D *h_ratio;
    TH1D *h_mc_ext_sum;

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c = new TCanvas();
    THStack *h_stack = new THStack();

    for (unsigned int k = 0; k < _util.type_prefix.size(); k++) {
        _util.GetHist(f_nuexsec, hist.at(k), Form("Flash/%s_%s", histname.c_str(), _util.type_prefix.at(k).c_str()));
        
        if (hist.at(k) == NULL) {
            std::cout << "Couldn't get all the flash histograms so exiting function..." << std::endl;
            return;
        }
    }

    topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    _util.SetTPadOptions(topPad, bottomPad);

    for (unsigned int i = 0; i < hist.size(); i++) {

        if (i == _util.k_data) {

            hist.at(i)->SetStats(kFALSE);
            hist.at(i)->SetMarkerStyle(20);
            hist.at(i)->SetMarkerSize(0.5);
        }

        // Scale EXT
        else if (i == _util.k_ext) {

            hist.at(i)->SetStats(kFALSE);
            // hist.at(i)->Scale(_util.ext_scale_factor);
        }

        // Scale Dirt
        else if (i == _util.k_dirt) {

            hist.at(i)->SetStats(kFALSE);
            // hist.at(i)->Scale(d_util.irt_scale_factor);
            hist.at(i)->SetFillColor(2);
            hist.at(i)->SetFillStyle(3354);
            hist_integrals.at(i) = hist.at(i)->Integral();
            integral_mc_ext += hist_integrals.at(_util.k_dirt);
        }

        // Scale MC
        else {
            hist.at(i)->SetStats(kFALSE);
            // hist.at(i)->Scale(_util.mc_scale_factor);
            hist.at(i)->SetFillColor(30);
            hist_integrals.at(i) = hist.at(i)->Integral();
            integral_mc_ext += hist_integrals.at(i);
        }
    }

    // Subtract the off beam from the on beam
    hist.at(_util.k_data)->Add(hist.at(_util.k_ext), -1);
    hist_integrals.at(_util.k_data) = hist.at(_util.k_data)->Integral();

    // Add the histograms to the stack
    h_stack->Add(hist.at(_util.k_mc));
    h_stack->Add(hist.at(_util.k_dirt));

    h_stack->Draw("hist");
    hist.at(_util.k_data)->Draw("same PE");

    // h_stack->GetYaxis()->SetTitle("Entries");

    // MC error histogram ------------------------------------------------------
    TH1D *h_error_hist = (TH1D *)hist.at(_util.k_mc)->Clone("h_error_hist");
    h_error_hist->Add(hist.at(_util.k_dirt), 1);

    h_error_hist->SetFillColorAlpha(12, 0.15);
    h_error_hist->Draw("e2, same");

    TLegend *leg_stack = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(hist.at(_util.k_data), "Beam-On - Beam-Off Data", "lep");
    leg_stack->AddEntry(hist.at(_util.k_dirt), "Out-of Cryo", "f");
    leg_stack->AddEntry(hist.at(_util.k_mc),   "Overlay", "f");

    leg_stack->Draw();

    bottomPad->cd();

    h_ratio = (TH1D *)hist.at(_util.k_data)->Clone("h_ratio");
    h_mc_ext_sum = (TH1D *)hist.at(_util.k_mc)->Clone("h_mc_ext_sum");

    // Add the dirt to overlay
    h_mc_ext_sum->Add(hist.at(_util.k_dirt), 1);

    h_ratio->Add(h_mc_ext_sum, -1);
    h_ratio->Divide(h_mc_ext_sum);

    h_ratio->GetXaxis()->SetLabelSize(12);
    h_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetYaxis()->SetLabelSize(11);
    h_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_ratio->GetXaxis()->SetTitleOffset(3.0);
    h_ratio->GetXaxis()->SetTitleSize(17);
    h_ratio->GetXaxis()->SetTitleFont(46);
    h_ratio->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);

    h_ratio->GetYaxis()->SetRangeUser(-0.5, 0.5);
    h_ratio->GetYaxis()->SetTitle("(Data - MC) / MC ");
    h_ratio->GetYaxis()->SetTitleSize(13);
    h_ratio->GetYaxis()->SetTitleFont(44);
    h_ratio->GetYaxis()->SetTitleOffset(1.5);
    h_ratio->SetTitle(" ");
    h_ratio->Draw("E");

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Draw Data to MC ratio
    _util.Draw_Data_MC_Ratio(c, double(hist_integrals.at(_util.k_data) * 1.0 / integral_mc_ext * 1.0), 0.34, 0.936, 0.34, 0.936);

    // Draw other data specifc quantities
    _util.Draw_Data_POT(c, Data_POT, 0.45, 0.915, 0.45, 0.915);

    // Add the weight labels
    Draw_WeightLabels(c);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeEfficiencyPlot(const char *print_name) {

    TTree *mc_tree;

    TFile *f_mc;

    // The uglyest hardcoded monstosity known to the universe. SORT this out KRISH...
    f_mc = TFile::Open(Form("files/trees/nuexsec_selected_tree_mc_run%s.root", _util.run_period));

    std::vector<double> efficiency_v; // efficiency vector
    std::vector<double> eff_err_v;    // efficiency error vector
    std::vector<double> purity_v;     // purity vector

    double efficiency, purity, eff_err;

    _util.GetTree(f_mc, mc_tree, "mc_eff_tree");
    mc_tree->SetBranchAddress("efficiency", &efficiency);
    mc_tree->SetBranchAddress("eff_err", &eff_err);
    mc_tree->SetBranchAddress("purity", &purity);

    int num_entries = mc_tree->GetEntries();

    // Fill the efficiency and purity vectors
    for (int y = 0; y < num_entries; y++) {
        mc_tree->GetEntry(y);
        efficiency_v.push_back(efficiency);
        eff_err_v.push_back(eff_err);
        purity_v.push_back(purity);
    }

    TCanvas *c = new TCanvas("c", "c", 600, 500);
    TH1D *h_eff = new TH1D("h_efficiency", "", efficiency_v.size(), 0, efficiency_v.size());
    TH1D *h_pur = new TH1D("h_purity", "", efficiency_v.size(), 0, efficiency_v.size());

    // c->SetGrid();
    c->SetGridy();
    // c->SetRightMargin(0.2);

    TLegend *leg_stack = new TLegend(0.7, 0.91, 0.90, 0.72);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    for (unsigned int k = 0; k < efficiency_v.size(); k++) {
        h_eff->Fill(_util.cut_dirs_pretty.at(k).c_str(), efficiency_v.at(k));
        h_pur->Fill(_util.cut_dirs_pretty.at(k).c_str(), purity_v.at(k));
        h_eff->SetBinError(k + 1, 0);
        h_pur->SetBinError(k + 1, 0);
    }
    
    TH1D* h_eff_clone = (TH1D*)h_eff->Clone();

    // Set the error bar on the efficiency histogram
    for (int bin = 0; bin < h_eff->GetNbinsX(); bin++) {
        std::cout <<eff_err_v.at(bin) << "  " << efficiency_v.at(bin) << std::endl;
        h_eff_clone->SetBinError(bin+1, eff_err_v.at(bin));
    }

    leg_stack->AddEntry(h_eff, "Efficiency", "ELP");
    leg_stack->AddEntry(h_pur, "Purity", "lp");

    h_eff->GetYaxis()->SetRangeUser(0.0, 1.1);
    h_eff->SetStats(kFALSE);
    h_eff->SetMarkerStyle(20);
    h_eff->SetMarkerSize(0.5);
    h_eff->SetLineWidth(2);
    h_eff->GetXaxis()->SetTitleSize(0.05);
    h_eff->GetYaxis()->SetLabelSize(0.05);
    h_eff->GetYaxis()->SetTitleSize(0.05);
    gPad->SetBottomMargin(0.12);
    h_eff->GetXaxis()->SetLabelFont(62);
    h_eff->GetXaxis()->SetLabelSize(0.03);

    h_eff->Draw("LP");

    h_pur->SetLineColor(kRed + 2);
    h_pur->SetStats(kFALSE);
    h_pur->SetMarkerStyle(20);
    h_pur->SetMarkerSize(0.5);
    h_pur->SetLineWidth(2);
    h_pur->Draw("LP,same");

    leg_stack->Draw();

    // Draw vertical lines to help the eye
    TLine *line;
    for (unsigned int l = 1; l < efficiency_v.size() + 1; l++) {
        line = new TLine(h_eff->GetBinCenter(l), 0, h_eff->GetBinCenter(l), 1.1);
        line->SetLineColor(12);
        line->SetLineStyle(kDotted);
        line->Draw();
    }

    h_eff_clone->Draw("E,X0,same");

    h_eff->GetXaxis()->SetTickLength(0.00);
    h_pur->GetXaxis()->SetTickLength(0.00);

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

    c->Print(print_name);

    delete c;
}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeEfficiencyPlotByCut(std::string var, bool mask_title) {

    std::vector<TH1D *> hist(_util.k_cuts_MAX); // The vector of histograms from the file for the plot
    std::vector<TEfficiency *> TEff_v(_util.k_cuts_MAX);

    TH1D *h_clone;

    // Helps determine what axes labels to draw 
    std::string var_string;

    if (var == "h_true_nu_E" || var == "h_true_nu_E_single_bin")              var_string = "nu";
    else if (var == "h_true_elec_E")       var_string = "elec";
    else if (var == "h_true_elec_E_rebin") var_string = "elec";
    else if (var == "h_true_nu_E_nue")     var_string = "nue";
    else if (var == "h_true_nu_E_nuebar")  var_string = "nuebar";
    else if (var == "h_true_nu_E_nue_single_bin")     var_string = "nue";
    else if (var == "h_true_nu_E_nuebar_single_bin")  var_string = "nuebar";
    else return;

    // Loop over the classifications and get the histograms
    for (unsigned int i = 0; i < _util.k_cuts_MAX; i++) {

        // MC
        _util.GetHist(f_nuexsec, hist.at(i), Form("TEff/%s_%s", var.c_str() ,_util.cut_dirs.at(i).c_str()));
        
        if (hist.at(i) == NULL)
            return;
    }

    // Loop over the cuts and draw the efficiencies
    for (int p = 0; p < _util.k_cuts_MAX; p++) {

        TCanvas * c = new TCanvas(Form("c_eff_by_cut_%s_%s_%s", _util.run_period,_util.cut_dirs.at(p).c_str(), var_string.c_str()), "c", 500, 500);
        c->SetTopMargin(0.11);

        h_clone = (TH1D *)hist.at(p)->Clone("h_clone");
        
        std::vector<double> eff_err;
        eff_err.resize(h_clone->GetNbinsX());
        
        // Get the bin errors based on binomial dist = sqrt(e/N*(1-e))) where e = n/N is the efficiency
        for (int bin = 0; bin < h_clone->GetNbinsX(); bin++){
            double n = h_clone->GetBinContent(bin+1); // selected = n
            double N = hist.at(_util.k_unselected)->GetBinContent(bin+1); // generated = N
            eff_err.at(bin) = (1.0/std::sqrt(N))*std::sqrt( (n/N)*(1.0 - (n/N)));
        }
        
        h_clone->Divide(hist.at(_util.k_unselected));
        
        // Now set the bin errors after the divide
        for (int bin = 0; bin < h_clone->GetNbinsX(); bin++){
            h_clone->SetBinError(bin+1, eff_err.at(bin));
        }
        
        
        h_clone->SetStats(kFALSE);
        
        if (var_string == "nu")     h_clone->SetTitle(Form("%s;True #nu_{e} + #bar{#nu}_{e} Energy [GeV]; Efficiency", _util.cut_dirs_pretty.at(p).c_str()));
        if (var_string == "elec")   h_clone->SetTitle(Form("%s;True Electron Energy [GeV]; Efficiency", _util.cut_dirs_pretty.at(p).c_str()));
        if (var_string == "nue")    h_clone->SetTitle(Form("%s;True #nu_{e} Energy [GeV]; Efficiency", _util.cut_dirs_pretty.at(p).c_str()));
        if (var_string == "nuebar") h_clone->SetTitle(Form("%s;True #bar{#nu}_{e} Energy [GeV]; Efficiency", _util.cut_dirs_pretty.at(p).c_str()));

        h_clone->GetYaxis()->SetRangeUser(0, 1);
        h_clone->SetLineColor(kBlack);
        h_clone->SetLineWidth(2);
        _util.IncreaseLabelSize(h_clone, c);
        if (mask_title) h_clone->SetTitle("");
        h_clone->Draw("E same");

        TH1D *h_true_nue = (TH1D *)hist.at(_util.k_unselected)->Clone("h_clone_true");

        gPad->SetRightMargin(0.17);

        Float_t rightmax = 1.1 * h_true_nue->GetMaximum();
        Float_t scale = gPad->GetUymax() / rightmax;
        h_true_nue->SetLineColor(kAzure - 6);
        h_true_nue->SetLineWidth(2);
        h_true_nue->Scale(scale);
        h_true_nue->Draw("hist,same");

        c->Update();

        TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax, 510, "+L");
        if (var_string == "nu")axis->SetTitle("True #nu_{e} + #bar{#nu}_{e} Events in FV");
        if (var_string == "elec")axis->SetTitle("True Electron Events in FV");
        if (var_string == "nue") axis->SetTitle("True #nu_{e} Events in FV");
        if (var_string == "nuebar") axis->SetTitle("True #bar{#nu}_{e} Events in FV");
        axis->SetTitleOffset(1.8);
        axis->SetLineColor(kAzure - 6);
        axis->SetLabelColor(kAzure - 6);
        axis->SetTitleColor(kAzure - 6);
        axis->SetTextFont(42);
        axis->SetLabelFont(42);
        axis->Draw();

        // Draw the run period on the plot
        _util.Draw_Run_Period(c, 0.76, 0.915, 0.76, 0.915);

        if (var == "h_true_nu_E")         c->Print(Form("plots/run%s/Efficiency/TEff_%s_nu_E.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));
        if (var == "h_true_elec_E")       c->Print(Form("plots/run%s/Efficiency/TEff_%s_elec_E.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));
        if (var == "h_true_elec_E_rebin") c->Print(Form("plots/run%s/Efficiency/TEff_%s_elec_E_rebin.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));
        if (var == "h_true_nu_E_single_bin") c->Print(Form("plots/run%s/Efficiency/TEff_%s_nu_E_single_bin.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));
        if (var == "h_true_nu_E_nue") c->Print(Form("plots/run%s/Efficiency/TEff_%s_nu_E_nue.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));
        if (var == "h_true_nu_E_nuebar") c->Print(Form("plots/run%s/Efficiency/TEff_%s_nu_E_nuebar.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));
        if (var == "h_true_nu_E_nue_single_bin") c->Print(Form("plots/run%s/Efficiency/TEff_%s_nu_E_nue_single_bin.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));
        if (var == "h_true_nu_E_nuebar_single_bin") c->Print(Form("plots/run%s/Efficiency/TEff_%s_nu_E_nuebar_single_bin.pdf", _util.run_period, _util.cut_dirs.at(p).c_str()));

        delete c;
    }
}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeInteractionPlot(const char *print_name, std::string cut_type, std::string flav, bool scale){

    std::vector<TH1D *> hist(_util.interaction_types.size());

    TCanvas * c = new TCanvas(Form("c_%s", print_name), "c", 500, 500);
    THStack *h_stack = new THStack();

    c->SetTopMargin(0.11);

    // Derive scale factor to scale the interaction plots to 6.0e20 POT
    double scale_factor{1.0};
    if (scale){
        if (std::string(_util.run_period) =="1") 
            scale_factor = 6.0e20 / _util.config_v.at(_util.k_Run1_MC_POT);
        else if (std::string(_util.run_period) == "3")
            scale_factor = 6.0e20 / _util.config_v.at(_util.k_Run3_MC_POT);
        else
            std::cout << "Unknown Run Period Configured" << std::endl;
    }


    for (unsigned int k = 0; k < hist.size(); k++) {
        
        _util.GetHist(f_nuexsec, hist.at(k), Form("Interaction/h_true_%s_E_%s_%s", flav.c_str(), _util.interaction_types.at(k).c_str(), cut_type.c_str()));
        

        if (hist.at(k) == NULL) {
            std::cout << "Couldn't get all the interaction histograms so exiting function..." << std::endl;
            return;
        }

        hist.at(k)->Scale(scale_factor);
    }

    hist.at(_util.k_plot_qe) ->SetFillColor(30);
    hist.at(_util.k_plot_res)->SetFillColor(38);
    hist.at(_util.k_plot_dis)->SetFillColor(28);
    hist.at(_util.k_plot_coh)->SetFillColor(42);
    hist.at(_util.k_plot_mec)->SetFillColor(36);
    hist.at(_util.k_plot_nc) ->SetFillColor(1);

    // Add the histograms to the stack
    for (unsigned int k = 0; k < hist.size(); k++) {
        h_stack->Add(hist.at(k));
    }

    h_stack->Draw("hist");
    h_stack->GetXaxis()->SetLabelSize(0.05);
    h_stack->GetXaxis()->SetTitleSize(0.05);
    h_stack->GetYaxis()->SetLabelSize(0.05);
    h_stack->GetYaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
    c->Update();
    if (flav == "nue_nuebar")h_stack->GetXaxis()->SetTitle("True #nu_{e} + #bar{#nu}_{e} Energy");
    if (flav == "nue")h_stack->GetXaxis()->SetTitle("True #nu_{e} Energy");
    if (flav == "nuebar")h_stack->GetXaxis()->SetTitle("True #bar{#nu}_{e} Energy");
    h_stack->GetYaxis()->SetTitle("Entries");

    if (scale) h_stack->SetMaximum(625);
    if (flav == "nuebar") h_stack->SetMaximum(100);

    TH1D *h_sum = (TH1D *)hist.at(_util.k_plot_qe)->Clone("h_sum");

    for (unsigned int i = 1; i < hist.size(); i++) {
        h_sum->Add(hist.at(i), 1);
    }

    h_sum->Draw("same E");

    TLegend *leg_stack = new TLegend(0.65, 0.89, 0.87, 0.6);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    // leg_stack->AddEntry(hist.at(_util.k_plot_nc), "NC", "f");
    leg_stack->AddEntry(hist.at(_util.k_plot_mec), "CC MEC", "f");
    leg_stack->AddEntry(hist.at(_util.k_plot_coh), "CC Coh", "f");
    leg_stack->AddEntry(hist.at(_util.k_plot_dis), "CC DIS", "f");
    leg_stack->AddEntry(hist.at(_util.k_plot_res), "CC Res", "f");
    leg_stack->AddEntry(hist.at(_util.k_plot_qe),  "CC QE", "f");

    leg_stack->Draw();

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
void HistogramPlotter::MakeInteractionEfficiency(const char *print_name){

    gStyle->SetOptStat(0);

    std::vector<std::vector<TH1D*>> hist(2);

    hist.at(0).resize(_util.interaction_types.size());  // Unselected
    hist.at(1).resize(_util.interaction_types.size());  // Selected

    TCanvas *c = new TCanvas();

    for (unsigned int k = 0; k < hist.at(0).size(); k++) {
        
        _util.GetHist(f_nuexsec, hist.at(0).at(k), Form("Interaction/h_true_nue_nuebar_E_%s_unselected", _util.interaction_types.at(k).c_str()));
        _util.GetHist(f_nuexsec, hist.at(1).at(k), Form("Interaction/h_true_nue_nuebar_E_%s_selected",   _util.interaction_types.at(k).c_str()));
        
        if (hist.at(0).at(k) == NULL || hist.at(1).at(k) == NULL) {
            std::cout << "Couldn't get all the interaction histograms so exiting function..." << std::endl;
            return;
        }

        hist.at(0).at(k) ->GetXaxis()->SetTitle("True #nu_{e} Energy");
        hist.at(0).at(k) ->GetYaxis()->SetTitle("Efficiency");

        hist.at(1).at(k) ->GetXaxis()->SetTitle("True #nu_{e} Energy");
        hist.at(1).at(k) ->GetYaxis()->SetTitle("Efficiency");

    }

    std::vector<TH1D*> h_ratio ;
    h_ratio.resize(_util.interaction_types.size());
    
    // Make the ratio histogram
    for (unsigned int type = 0; type < hist.at(0).size(); type++){
        
        h_ratio.at(type) = (TH1D*) hist.at(1).at(type)->Clone(Form("h_ratio_%s" ,_util.interaction_types.at(type).c_str()) );
        h_ratio.at(type)  ->Divide(hist.at(0).at(type));

        if (type == _util.k_plot_qe){
            h_ratio.at(_util.k_plot_qe) ->SetLineColor(30); 
        }
        if (type == _util.k_plot_res){
            h_ratio.at(_util.k_plot_res)->SetLineColor(38);
        }
        if (type == _util.k_plot_dis){
            h_ratio.at(_util.k_plot_dis)->SetLineColor(28);
        }
        if (type == _util.k_plot_coh){
            h_ratio.at(_util.k_plot_coh)->SetLineColor(42);
        }
        if (type == _util.k_plot_mec){
            h_ratio.at(_util.k_plot_mec)->SetLineColor(36);
        }


        _util.IncreaseLabelSize(h_ratio.at(type), c);
        h_ratio.at(type) ->SetFillColor(0);
        h_ratio.at(type)->GetYaxis()->SetRangeUser(0, 0.6);
        h_ratio.at(type)->GetXaxis()->SetRangeUser(0, 3.5);
        h_ratio.at(type)->SetLineWidth(2);

        if (type == _util.k_plot_nc || type == _util.k_plot_coh) continue; // Too low stats for the plot
        h_ratio.at(type)->Draw("hist,L,same");

    }

    
    

    TLegend *leg_stack = new TLegend(0.55, 0.89, 0.75, 0.6);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    // leg_stack->AddEntry(hist.at(_util.k_plot_nc), "NC", "f");
    leg_stack->AddEntry(h_ratio.at(_util.k_plot_mec), "CC MEC", "l");
    // leg_stack->AddEntry(h_ratio.at(_util.k_plot_coh), "CC Coh", "l"); // Too low stats for the plot
    leg_stack->AddEntry(h_ratio.at(_util.k_plot_dis), "CC DIS", "l");
    leg_stack->AddEntry(h_ratio.at(_util.k_plot_res), "CC Res", "l");
    leg_stack->AddEntry(h_ratio.at(_util.k_plot_qe),  "CC QE",  "l");

    leg_stack->Draw();

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
void HistogramPlotter::Plot2D_Signal_Background(const char *print_name, const char *histname){

    std::vector<TH2D *> hist(_util.sig_bkg_prefix.size());

    TCanvas * c = new TCanvas(Form("c_%s", print_name), "c", 500, 500);
    c->SetTopMargin(0.11);

    for (unsigned int k = 0; k < hist.size(); k++){
        
        _util.GetHist(f_nuexsec, hist.at(k), Form("2D/%s_%s", histname, _util.sig_bkg_prefix.at(k).c_str()));
        
        if (hist.at(k) == NULL){
            std::cout << "Couldn't get all the interaction histograms so exiting function..." << std::endl;
            return;
        }
        
        hist.at(k)->SetStats(kFALSE);
    }

    TLegend *leg_stack = new TLegend(0.6, 0.89, 0.87, 0.75);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(hist.at(_util.k_signal), "Signal", "f");
    leg_stack->AddEntry(hist.at(_util.k_background), "Background", "f");

    hist.at(_util.k_background)->Draw("box");
    hist.at(_util.k_signal)->Draw("box, same");
    // hist.at(_util.k_background) ->SetFillStyle(3244);
    hist.at(_util.k_background)->SetFillColorAlpha(kRed + 2, 0.2);
    hist.at(_util.k_background)->Draw("box,same");

    // IncreaseLabelSize(hist.at(_util.k_background));

    // Draw cut lines to help the eye
    std::vector<TLine *> line_v;
    line_v.resize(9);
    line_v.at(0) = new TLine(0, 3, 1.75, 3);

    line_v.at(1) = new TLine(1.75, 3, 1.75, 12);
    line_v.at(2) = new TLine(1.75, 12, 2.5, 12);

    line_v.at(3) = new TLine(2.5, 3, 2.5, 12);
    line_v.at(4) = new TLine(2.5, 3, 3.5, 3);

    line_v.at(5) = new TLine(3.5, 0, 3.5, 3);
    line_v.at(6) = new TLine(3.5, 0, 4.5, 0);

    line_v.at(7) = new TLine(4.5, 0, 4.5, 3);
    line_v.at(8) = new TLine(4.5, 3, 10.0, 3);

    for (unsigned int l = 0; l < line_v.size(); l++) {
        line_v.at(l)->SetLineColor(kBlack);
        line_v.at(l)->SetLineWidth(2);
        line_v.at(l)->Draw();
    }

    leg_stack->Draw();

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
void HistogramPlotter::Save1DHists(const char *print_name, const char *histname, std::string cut_type, bool scale) {

    TH1D *hist;
    _util.GetHist(f_nuexsec, hist, Form("True/%s_MC_%s", histname, cut_type.c_str()));

    if (hist == NULL)
        std::cout << "couldn't get the hist!" << std::endl;

    // Derive scale factor to scale the interaction plots to 6.0e20 POT
    double scale_factor{1.0};
    if (scale){
        if (std::string(_util.run_period) =="1") 
            scale_factor = 6.0e20 / _util.config_v.at(_util.k_Run1_MC_POT);
        else if (std::string(_util.run_period) == "3")
            scale_factor = 6.0e20 / _util.config_v.at(_util.k_Run3_MC_POT);
        else
            std::cout << "Unknown Run Period Configured" << std::endl;
    }

    if (scale) hist->Scale(scale_factor);

    TCanvas * c = new TCanvas(Form("c_%s", print_name), "c", 500, 500);
    c->SetTopMargin(0.11);

    hist->SetStats(kFALSE);

    _util.IncreaseLabelSize(hist, c);

    hist->SetLineColor(kAzure - 6);
    hist->SetLineWidth(2);
    hist->Draw("hist,E");

    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
void HistogramPlotter::Save2DHists(const char *print_name, const char *histname, std::string cut_type, bool yex) {

    TH2D *hist;
    _util.GetHist(f_nuexsec, hist, Form("True/%s_MC_%s", histname, cut_type.c_str()));

    if (hist == NULL)
        std::cout << "couldn't get the hist!" << std::endl;

    TCanvas * c = new TCanvas(Form("c_%s", print_name), "c", 500, 500);
    c->SetTopMargin(0.11);

    hist->SetStats(kFALSE);

    _util.IncreaseLabelSize(hist, c);

    hist->Draw("colz");

    c->Update();

    // If yex (y=x) draw a y=x line to guide the eye
    TLine * line;
    if (yex){
    
        line = new TLine(0, 0, gPad->GetUymax(), gPad->GetUymax());
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw();
    }
    
    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.76, 0.915, 0.76, 0.915);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
void HistogramPlotter::Save2DHistsNorm(const char *print_name, const char *histname, std::string cut_type, bool yex, std::string normtype) {

    TH2D *hist, *hist_clone;
    _util.GetHist(f_nuexsec, hist, Form("True/%s_MC_%s", histname, cut_type.c_str()));

    if (hist == NULL)
        std::cout << "couldn't get the hist!" << std::endl;

    hist_clone = (TH2D*)hist->Clone(Form("h_%s_clone", print_name));

    TCanvas * c = new TCanvas(Form("c_%s", print_name), "c", 500, 500);
    c->SetTopMargin(0.11);

    hist->SetStats(kFALSE);

    _util.IncreaseLabelSize(hist, c);


    // Now we normalise by column (true) or row (reco)
    
    // We are normalising by truth integral
    if (normtype == "true"){

        // Loop over rows
        for (int row=1; row<hist_clone->GetXaxis()->GetNbins()+1; row++) {
            double integral = 0;

            // Loop over columns and get the integral
            for (int col=1; col<hist_clone->GetYaxis()->GetNbins()+1; col++){
                integral+=hist_clone->GetBinContent(row, col);            
            }

            // Now normalise the column entries by the integral
            for (int col=1; col<hist_clone->GetYaxis()->GetNbins()+1; col++){
                hist_clone->SetBinContent(row,col, hist_clone->GetBinContent(row, col)/ integral );
                
            }
        } 
    }
    // We normalise by the row intgral (i.e. reco space)
    else {
        
        // Loop over rows
        for (int col=1; col<hist_clone->GetYaxis()->GetNbins()+1; col++) {
            double integral = 0;

            // Loop over columns and get the integral
            for (int row=1; row<hist_clone->GetXaxis()->GetNbins()+1; row++){
                integral+=hist_clone->GetBinContent(row, col);            
            }

            // Now normalise the column entries by the integral
            for (int row=1; row<hist_clone->GetXaxis()->GetNbins()+1; row++){
                hist_clone->SetBinContent(row,col, hist_clone->GetBinContent(row, col)/ integral );
                
            }
        } 
    }

    hist_clone->SetMinimum(0.0);
    hist_clone->SetMaximum(1.0);
    gStyle->SetPaintTextFormat("4.2f");
    hist_clone->SetMarkerSize(0.8);
    hist_clone->SetMarkerColor(kBlack);
    hist_clone->Draw("colz");

    TH2D* hist_clone2 = (TH2D*)hist_clone->Clone(Form("h_%s_clone2", print_name));

    TLatex* range;
    if (normtype == "true")range = new TLatex(0.5,0.92, "Column Normalised");
    else range = new TLatex(0.45,0.92,"Row Normalised");
    range->SetTextColor(kGray+2);
    range->SetNDC();
    range->SetTextSize(0.038);
    range->SetTextAlign(32);
    range->Draw();

    c->Update(); 

    // If yex (y=x) draw a y=x line to guide the eye
    TLine * line;
    if (yex){
    
        line = new TLine(0, 0, gPad->GetUymax(), gPad->GetUymax());
        line->SetLineColor(kGray);
        line->SetLineWidth(2);
        line->Draw();
    }

    hist_clone2->Draw("text,same");
    
    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.76, 0.915, 0.76, 0.915);

    c->Print(print_name);
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
