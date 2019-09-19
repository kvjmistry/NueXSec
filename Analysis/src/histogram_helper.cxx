#include "../include/histogram_helper.h"

// -----------------------------------------------------------------------------
histogram_helper::~histogram_helper() { 
    
    // Make sure the file is closed
    f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void histogram_helper::MakeDirectory(std::string type){
        
    f_nuexsec->cd();

    TDirectory *top_dir; // e.g MC, Data, EXT
    bool bool_dir;       // Check if directory exists already
    TString type_tstr = type;
   
    // Create the top directory
    bool_dir = _utility_instance.GetDirectory(f_nuexsec, top_dir, type_tstr );
    if (!bool_dir) top_dir = f_nuexsec->mkdir(type.c_str());
    
    // Make the the top dir the current directory
    top_dir->cd();

    // Create subdirectory for cut type
    TDirectory *dir_plot_types[plot_types.size()];
    
    // Create a new subdirectory for each cut
    const Int_t ncuts = cut_dirs.size();
    TDirectory *dir_cut[ncuts];

    // Create a new subdirectory for each classification
    TDirectory *dir_classification[classification_dirs.size()];
    
    // Loop over the plot types ------------------------------------------------
    for (unsigned int k = 0; k < plot_types.size(); k++) {
        
        // Get the directory 
        bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_plot_types[k] ,Form("%s/%s", type.c_str(), plot_types.at(k).c_str()) );

        // Make the directory
        if (!bool_dir) dir_plot_types[k] = top_dir->mkdir(plot_types.at(k).c_str());

        dir_plot_types[k]->cd();

        // If we have stacked histograms, we make plots by cut 
        if (plot_types.at(k) == "Stack"){
            
            // Loop over the cuts ----------------------------------------------
            for (int i = 0; i < ncuts; i++) {
               
                // Get the directory 
                bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_cut[i] ,Form("%s/%s/%s", type.c_str(), plot_types.at(k).c_str(), cut_dirs.at(i).c_str()));

                // Make the directory
                if (!bool_dir) dir_cut[i] = dir_plot_types[k]->mkdir(cut_dirs.at(i).c_str());
                dir_cut[i]->cd();
                
                // Loop over the classifications -------------------------------
                for (unsigned int j = 0; j < classification_dirs.size(); j++){
                
                    // Get the directory 
                    bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_classification[j] ,Form("%s/%s/%s/%s", type.c_str(), plot_types.at(k).c_str(), cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()));

                    // Make the directory
                    if (!bool_dir) dir_classification[j] = dir_cut[i]->mkdir(classification_dirs.at(j).c_str());
                    dir_classification[j]->cd();
                    top_dir->cd();    // change current directory to top

                } // End loop over the classifications -------------------------
                
                top_dir->cd();    // change current directory to top
                
            } // End loop over the cuts ----------------------------------------
       
        }
         
        top_dir->cd();    // change current directory to top
    
    } // End loop over plot types ----------------------------------------------

    top_dir->Write("",TObject::kOverwrite);
   
}
// -----------------------------------------------------------------------------
void histogram_helper::Initialise(){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject("nuexsec.root") ) {
        f_nuexsec = new TFile("nuexsec.root", "UPDATE");
    }

    MakeDirectory("MC");
    MakeDirectory("Data");
    MakeDirectory("EXT");
    MakeDirectory("Dirt");
}
// -----------------------------------------------------------------------------
void histogram_helper::InitHistograms(){
    
    // Flash Histograms
    h_flash_time_v.resize(k_flash_MAX);
    for (unsigned int i=0; i < h_flash_time_v.size();i++){
        h_flash_time_v.at(i) = new TH1D ( Form("h_flash_time_%s", type_prefix.at(i).c_str()) ,"", 80, 0, 20);
    }

    // Reco Vtx X
    h_reco_vtx_x.resize(k_cuts_MAX);
    h_reco_vtx_y.resize(k_cuts_MAX);
    h_reco_vtx_z.resize(k_cuts_MAX);
    
    for (unsigned int i=0; i < cut_dirs.size();i++){

        h_reco_vtx_x.at(i).resize(k_classifications_MAX);
        h_reco_vtx_y.at(i).resize(k_classifications_MAX);
        h_reco_vtx_z.at(i).resize(k_classifications_MAX);

        for (unsigned int j=0; j < classification_dirs.size();j++){
            h_reco_vtx_x.at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 270);
            h_reco_vtx_y.at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 120);
            h_reco_vtx_z.at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 1050);
        }
    }

}
// -----------------------------------------------------------------------------
void histogram_helper::FillMCTruth( double mc_nu_energy,  double mc_nu_momentum,  int mc_nu_id, bool in_tpc,
                                    double mc_nu_vtx_x,   double mc_nu_vtx_y,  double mc_nu_vtx_z,
                                    double mc_nu_dir_x,   double mc_nu_dir_y,  double mc_nu_dir_z,
                                    double mc_ele_dir_x,  double mc_ele_dir_y, double mc_ele_dir_z,
                                    double mc_ele_energy, double mc_ele_momentum ) {
    
    double mc_cos_theta     = -999;
    double theta            = -999;
    double mc_ele_cos_theta = -999;
    double mc_ele_theta     = -999;

    // Caclulate theta and cos(theta)
    if (mc_nu_momentum != 0) {
        mc_cos_theta = mc_nu_dir_z;
        theta        = acos(mc_cos_theta) * (180 / 3.1415);
    }
    
    // Caclulate theta and cos(theta) for the electron
    if (mc_ele_momentum != 0) {
        mc_ele_cos_theta = mc_ele_dir_z;
        mc_ele_theta     = acos(mc_ele_cos_theta) * (180/3.1415);
    }
    
    // Calculate Phi
    double phi        = atan2(mc_nu_dir_y, mc_nu_dir_x) * (180/3.1415);
    double mc_ele_phi = atan2(mc_ele_dir_y, mc_ele_dir_x) * (180/3.1415);
   
    if ((mc_nu_id == 1 || mc_nu_id == 5) && in_tpc == true){
        
        // 1D Hists
        h_nue_true_theta->Fill(theta);
        h_nue_true_phi  ->Fill(phi);
        
        // 2D Hists
        h_nue_true_theta_phi   ->Fill(phi, theta );
        h_nue_true_energy_theta->Fill(mc_nu_momentum, theta);
        h_nue_true_energy_phi  ->Fill(mc_nu_momentum, phi);
        
        h_ele_true_energy_theta->Fill(mc_ele_energy, mc_ele_theta);
        h_ele_true_energy_phi  ->Fill(mc_ele_energy, mc_ele_phi);
    }
    

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteMCTruth(std::string type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth

    // Get the truth directory and cd
    bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s", type.c_str(), "Truth") );
    if (bool_dir) truth_dir->cd();

    // Now write the histograms
    h_nue_true_theta->Write("",TObject::kOverwrite);
    h_nue_true_phi  ->Write("",TObject::kOverwrite);
    
    // 2D Hists
    h_nue_true_theta_phi   ->Write("",TObject::kOverwrite);
    h_nue_true_energy_theta->Write("",TObject::kOverwrite);
    h_nue_true_energy_phi  ->Write("",TObject::kOverwrite);
    
    h_ele_true_energy_theta->Write("",TObject::kOverwrite);
    h_ele_true_energy_phi  ->Write("",TObject::kOverwrite);


}
// -----------------------------------------------------------------------------
void histogram_helper::FillOptical(std::vector<std::vector<double>> optical_list_flash_time_v, int type){

    // Loop over the optical list events
    for (unsigned int i = 0; i < optical_list_flash_time_v.size();i++){

        // Loop over the flashes in each event
        for (unsigned int j = 0; j < optical_list_flash_time_v.at(i).size();j++){
            h_flash_time_v.at(type)->Fill(optical_list_flash_time_v.at(i).at(j));
        }
    }

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteOptical(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth

    // Get the Optical directory and cd
    bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s", type_prefix.at(type).c_str(), "Optical") );
    if (bool_dir) truth_dir->cd();

    // Now write the histograms
    h_flash_time_v.at(type)->Write("",TObject::kOverwrite);
    
   
}
// -----------------------------------------------------------------------------
int histogram_helper::IndexOfClassification(std::string tpco_id){

    // Nue CC
    if (tpco_id == "nue_cc_qe"  || tpco_id == "nue_bar_cc_qe"  ||
        tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res" ||
        tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis" || 
        tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh" || 
        tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec") {
        return k_nue_cc;
    }

    // Nue CC OOFV
    if (tpco_id == "nue_cc_out_fv"){
        return k_nue_cc_out_fv;
    }

    // NuMu CC
    if (tpco_id == "numu_cc_qe" || tpco_id == "numu_cc_res" ||
       tpco_id == "numu_cc_dis" || tpco_id == "numu_cc_coh" || 
       tpco_id == "numu_cc_mec" || tpco_id == "numu_cc_mixed"){
        return k_numu_cc;
    }

    // NC
    if (tpco_id == "nc"){
        return k_nc;
    }
    
    // NC pi0
    if (tpco_id == "nc_pi0"){
        return k_nc_pi0;
    }

    // Nue CC Mixed
    if (tpco_id == "nue_cc_mixed"){
        return k_nue_cc_mixed;
    }

    // Cosmic
    if (tpco_id == "cosmic"){
        return k_cosmic;
    }

    // Other Mixed
    if (tpco_id == "other_mixed"){
        return k_nc_mixed;
    }

    // Unmatched
    if (tpco_id == "unmatched" || tpco_id == "bad_reco"){
        return k_unmatched;
    }

    // Data
    if (tpco_id == "Data"){
        return  k_leg_data;
    }

    // EXT/ In time cosmics
    if (tpco_id == "EXT"){
        return  k_leg_ext;
    }

    // Dirt
    if (tpco_id == "Dirt"){
        return  k_leg_dirt;
    }

    std::cout << "Reached end of index of classifcation, this is BAD!!!! " << tpco_id << std::endl;
    return k_classifications_MAX;
}
// -----------------------------------------------------------------------------
void histogram_helper::FillRecoVtx(int classification_index, int cut_index, const xsecAna::TPCObjectContainer &tpc_obj){

    h_reco_vtx_x.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxX());
    h_reco_vtx_y.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxY());
    h_reco_vtx_z.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxZ());

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteRecoVtx(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth

    // loop over the cut directories
    for (unsigned int i = 0; i < cut_dirs.size(); i++){
        
        // loop over the classification directories
        for (unsigned int j = 0; j < classification_dirs.size(); j++){

            // Get the classification directory and cd
            bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s/%s/%s", type_prefix.at(type).c_str(), "Stack", cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str() ) );
            if (bool_dir) truth_dir->cd();

            // Now write the histogram
            h_reco_vtx_x.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_y.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_z.at(i).at(j)->Write("",TObject::kOverwrite);
        }

    }
    
}
// -----------------------------------------------------------------------------
void histogram_helper::SetStack(std::string hist_name, std::string cut_name, bool area_norm,  bool logy, const char* x_axis_name,
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
        if (i == k_leg_data) continue; // Dont use the data
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
std::vector <double> histogram_helper::Chi2Calc(TH1D * h_mc_ext, TH1D * h_data, const bool area_norm, const double return_norm){
    const int n_bins = h_mc_ext->GetNbinsX();

    const double f_1 = h_mc_ext->Integral();
    const double f_2 = h_data->Integral();

    //area normalised?
    TH1D * h_mc_ext_clone = (TH1D*)h_mc_ext->Clone("h_mc_ext_clone");
    TH1D * h_data_clone = (TH1D*)h_data->Clone("h_data_clone");
    if(!area_norm) {h_mc_ext_clone->Scale(f_2/f_1); }
    
    if(area_norm) {
        //this keeps them area normalised,
        //but at the original values, not 0->1
        //which messes with the chi2 and p calc
        h_mc_ext_clone->Scale(return_norm);
        h_data_clone->Scale(return_norm);

        const double f_1_adj = h_mc_ext->Integral();
        const double f_2_adj = h_data->Integral();
        h_mc_ext_clone->Scale(f_2_adj/f_1_adj);

    }
    //h_data_clone->Scale(1./f_2);

    std::vector <double> chi2;
    double chi2_val = 0;
    double n_mc_ext_val = 0;
    double n_data_val = 0;
    for( int i = 1; i < n_bins; i++) {
        const double n_mc_ext = h_mc_ext_clone->GetBinContent(i);
        const double n_data   = h_data_clone  ->GetBinContent(i);

        //don't calculate chi2 for bins where no comparison possible
        if(n_data == 0 || n_mc_ext == 0) { continue; }

        //chi2_val += (pow((n_mc_ext - n_data),2) / n_mc_ext);
        //chi2_val += (pow((n_data - n_mc_ext),2)) / (((n_data * f_2) / pow(f_2, 2)) + ((n_mc_ext * f_1) / pow(f_1, 2)));
        chi2_val += 2 * (n_mc_ext - n_data + (n_data * TMath::Log(n_data/n_mc_ext)));

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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------