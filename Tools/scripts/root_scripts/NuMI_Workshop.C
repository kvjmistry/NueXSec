/*
NuMI_Workshop.C: A root macro that is a simplified version of the SystematicsHelper.cxx function in 
                 Krish's NueXsec Module: https://github.com/kvjmistry/NueXSec

                 This script has been made to make systematic plots from a file that contains
                 weighted histograms intended for the NuMI Workshop Jan 2021.

                 To run do:

                 root -l -b -q NuMI_Workshop.C 

                 This will make some plots in a folder called "plots"
                 take a look at them!

*/

// Declare functions for macro -- see comments in the implementation of the functions
void InitialseCV(std::vector<std::vector<TH1D*>> &cv_hist_vec, TFile* f);
void SetLabelName(std::string label, std::string &label_up, std::string &label_dn);
void PlotUnisim(std::string sys_label, int var, std::string label_pretty, TFile* f);
void PlotMultisim(std::string label, int var, std::string label_pretty, int universes, TFile *f);
void NuMI_Workshop();
void PlotDetVar(std::string label, int var, int detvar_index, std::string label_pretty, TFile *f);

// --

// Define some enums to make selecting the variable a little easier
enum TH1D_xsec_var_vars {
    k_integrated,     // Integrated X-Section
    k_reco_el_E,      // Reconstructed electron energy
    k_true_el_E,      // True electron energy
    k_var_MAX
};

// enum for histogram types
enum TH1D_xsec_hist_vars {
    k_xsec_sel,     // Selected event histogram binned in energy
    k_xsec_bkg,     // Bkg event histogram binned in energy
    k_xsec_gen,     // Gen event histogram binned in energy
    k_xsec_sig,     // Sig event histogram binned in energy
    k_xsec_eff,     // Efficiency histogram binned in energy
    k_xsec_ext,     // EXT event histogram binned in energy
    k_xsec_dirt,    // Dirt event histogram binned in energy
    k_xsec_data,    // Data event histogram binned in energy
    k_xsec_mcxsec,  // MC Cross Section
    k_xsec_dataxsec,// Data Cross Section
    k_TH1D_xsec_MAX
};

// enum for up and dn
enum updn {k_up, k_dn};

// Define the variable vector 
std::vector<std::string> vars = {"integrated", "reco_el_E", "true_el_E" };

// Define cross section types histogram
std::vector<std::string> xsec_types = {"sel", "bkg", "gen", "sig", "eff", "ext", "dirt", "data", "mc_xsec", "data_xsec"};
std::vector<std::string> xsec_types_pretty = {"Selected", "Background", "Generated Signal", "Signal", "Efficiency", "Beam-Off", "Dirt", "Beam-On", "MC", "Data"};

std::vector<std::string> var_labels_xsec = {";;#nu_{e} + #bar{#nu}_{e} CC Cross Section [10^{-39} cm^{2}/nucleon]",
                                           ";Reco. Leading Shower Energy [GeV];#frac{d#sigma}{dE^{reco}_{e#lower[-0.5]{-} + e^{+}}} [10^{-39} cm^{2}/GeV/nucleon]",
                                           ";True e#lower[-0.5]{-} + e^{+} Energy [GeV];#frac{d#sigma}{dE^{true}_{e#lower[-0.5]{-} + e^{+}}} [10^{-39} cm^{2}/GeV/nucleon"
                                           };

std::vector<std::string> var_labels_events = {";;Entries",
                                    ";Reco. Leading Shower Energy [GeV]; Entries / GeV",
                                    ";True e#lower[-0.5]{-} + e^{+} Energy [GeV]; Entries / GeV"
                                    };

std::vector<std::string> var_labels_eff = {";;Efficiency",
                                    ";Reco. Leading Shower Energy [GeV]; Efficiency",
                                    ";True e#lower[-0.5]{-} + e^{+} Energy [GeV]; Efficiency"
                                    };

// Containter for the central value histograms
std::vector<std::vector<TH1D*>> cv_hist_vec; // reco elec e, <gen, sig, etc>


// variation names for detvar
std::vector<std::string> detvar_string = {
        "CV",
        "LYRayleigh",
        "LYDown",
        "SCE",
        "Recomb2",
        "WireModX",
        "WireModYZ",
        "WireModThetaXZ",
        "WireModThetaYZ_withSigmaSplines"
    };

// strings used for the legend of the plots of the detector variations
std::vector<std::string> detvar_string_pretty = {
    "CV",
    "LY Rayleigh",
    "LY Down",
    "SCE",
    "Recombination",
    "WM X",
    "WM YZ",
    "WM Theta XZ",
    "WM Theta YZ w/ Spl."
};

enum enum_variations {
    k_CV,
    k_LYRayleigh,
    k_LYDown,
    k_SCE,
    k_Recomb2,
    k_WireModX,
    k_WireModYZ,
    k_WireModThetaXZ,
    k_WireModThetaYZ_withSigmaSplines,
    k_vars_MAX
};

// POT of the detvar
std::vector<double> POT_v = {
 3.68983e+20, // POT CV
 3.68926e+20, // POT LYRayleigh
 3.566e+20,   // POT LYDown
 3.71139e+20, // POT SCE
 3.72109e+20, // POT Recomb2
 3.6749e+20,  // POT WireModX
 3.69849e+20, // POT WireModYZ
 3.68008e+20, // POT WireModThetaXZ
 3.68995e+20  // POT WireModThetaYZ_withSigmaSplines
};

// ----------------------------      MAIN     -----------------------------------
// Main function that gets exectuted
void NuMI_Workshop(){

    // Who cares about stats?
    gStyle->SetOptStat(0);

    // Make the plots directory
    gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots; fi");

    // Get the input root histogram
    TFile *file_in = TFile::Open("crosssec_run1.root", "READ");

    // Grab the CV histograms from the file
    InitialseCV(cv_hist_vec, file_in);

    // Unisim Plot
    // Next try and make another unisim plot
    PlotUnisim("RPA", k_reco_el_E, "RPA", file_in); 

    // Multisim Plot
    // Next try and make the PPFX systematics plot or the Geant Reinteractions
    PlotMultisim("weightsGenie", k_reco_el_E,  "GENIE", 500, file_in);   

    // Detvar Unisim Plot
    // Next try and make another detvar unisim plot
    PlotDetVar("LYRayleigh", k_reco_el_E, k_LYRayleigh, detvar_string_pretty.at(k_LYRayleigh), file_in);

}
// -----------------------------------------------------------------------------
// Get the CV histograms from the file to plot against
void InitialseCV(std::vector<std::vector<TH1D*>> &cv_hist_vec, TFile* f){
    
    // Resize the vector of histograms
    cv_hist_vec.resize(vars.size());
    
    for (unsigned int var = 0; var < vars.size(); var++){
        cv_hist_vec.at(var).resize(xsec_types.size());
    }


    // Loop over the vars
    for (unsigned int var = 0; var < vars.size(); var++){
        
        // Loop over the typrs
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            cv_hist_vec.at(var).at(k) = (TH1D*)f->Get(Form( "CV/%s/h_run1_CV_0_%s_%s", vars.at(var).c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

            if (cv_hist_vec.at(var).at(k) == NULL)  
                std::cout << "Failed to get the histogram!" << std::endl;

            // Customise
            cv_hist_vec.at(var).at(k)->SetLineWidth(2);
            cv_hist_vec.at(var).at(k)->SetLineColor(kBlack);

            // Set the Titles
            if (k == k_xsec_mcxsec)
                cv_hist_vec.at(var).at(k)->SetTitle(var_labels_xsec.at(var).c_str());
            else if (k == k_xsec_dataxsec)
                cv_hist_vec.at(var).at(k)->SetTitle(var_labels_events.at(var).c_str());
            else if (k == k_xsec_eff)
                cv_hist_vec.at(var).at(k)->SetTitle(var_labels_eff.at(var).c_str());
            else
                cv_hist_vec.at(var).at(k)->SetTitle(var_labels_events.at(var).c_str());

        }

    }
}
// -----------------------------------------------------------------------------
// Function to handle the up/down labels in the file
void SetLabelName(std::string label, std::string &label_up, std::string &label_dn){

    if  (label == "Horn_curr"          ){
        label_up = "Horn_p2kA";
        label_dn = "Horn_m2kA";
    }
    else if  (label == "Horn1_x"       ){
        label_up = label + "_p3mm";
        label_dn = label + "_m3mm";
    }
    else if  (label == "Horn1_y"       ){
        label_up = label + "_p3mm";
        label_dn = label + "_m3mm";
    }
    else if  (label == "Beam_spot"    ){
        label_up = label + "_1_1mm";
        label_dn = label + "_1_5mm";
    }
    else if  (label == "Horn2_x"       ){
        label_up = label + "_p3mm";
        label_dn = label + "_m3mm";
    }
    else if  (label == "Horn2_y"       ){
        label_up = label + "_p3mm";
        label_dn = label + "_m3mm";
    }
    else if  (label == "Horn_Water"    ){
        label_up = "Horns_0mm_water";
        label_dn = "Horns_2mm_water";
    }
    else if  (label == "Beam_shift_x"  ){
        label_up = label + "_p1mm";
        label_dn = label + "_m1mm";
    }
    else if  (label == "Beam_shift_y"  ){
        label_up = label + "_p1mm";
        label_dn = label + "_m1mm";
    }
    else if  (label == "Target_z"      ){
        label_up = label + "_p7mm";
        label_dn = label + "_m7mm";
    }
    else if  (label == "Horn1_refined_descr"                || label == "Decay_pipe_Bfield" || label == "Old_Horn_Geometry" ||
              label == "LYRayleigh"                         || label == "LYDown"            || label == "LYAttenuation"     || label == "SCE" || label == "Recomb2"       || 
              label == "WireModX"                           || label == "WireModYZ"         || label == "WireModThetaXZ"    || label == "WireModThetaYZ_withSigmaSplines" || 
              label == "WireModThetaYZ_withoutSigmaSplines" || label == "WireModdEdX"       || label == "pi0"){
        label_up = label;
        label_dn = label;
    }
    // Other variation
    else {
        label_up = label + "up";
        label_dn = label + "dn";
    }

}
// -----------------------------------------------------------------------------
/*
    // Function to make a unisim plot:

    std::string sys_label - The name of the systematic in the file. Choices are: 
    
    RPA, CCMEC, AxFFCCQE, VecFFCCQE, DecayAngMEC, ThetaDelta2Npi, ThetaDelta2NRad,
    NormCCCOH, NormNCCOH, Dirt, POT, pi0, Horn1_x, Horn_curr, Horn1_y, Beam_spot, 
    Horn2_x, Horn2_y, Horn_Water, Beam_shift_x, Beam_shift_y, Target_z
    ---
    int var - the variable to plot the distributions as a function of. Choices are:
    
    k_integrated, k_reco_el_E, k_true_el_E
    ---
    std::string label_pretty -- name that goes in the legend, choose something appropriate:
    e.g if plotting RPA, then use RPA or Random Phase Approximation
    ---
    TFile* f -- the input TFile
*/
void PlotUnisim(std::string sys_label, int var, std::string label_pretty, TFile* f){

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c;

    std::vector<std::vector<TH1D*>> h_universe;
    
    // Resize to the number of universes
    h_universe.resize(2);
    h_universe.at(k_up).resize(xsec_types.size());
    h_universe.at(k_dn).resize(xsec_types.size());

    // Now get the histograms
    std::string label_up = sys_label + "up";
    std::string label_dn = sys_label + "dn";

    // Set the Unisim up down variation name
    SetLabelName(sys_label, label_up, label_dn);

    // Check if its just a single on/off type variation
    // This case we dont want to plot the up/dn, but just once
    bool single_var = false;
    if (label_up == label_dn) single_var = true;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        

        // Get the up variation ----
        h_universe.at(k_up).at(k) = (TH1D*)f->Get(Form( "%s/%s/h_run1_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        // Get the dn variation ----

        h_universe.at(k_dn).at(k) = (TH1D*)f->Get(Form( "%s/%s/h_run1_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
    }


    // Now we want to draw them
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        c         = new TCanvas("c", "c", 500, 500);
        topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        
        // Set options for the TPad and draw them
        topPad   ->SetBottomMargin(0.05);
        topPad   ->SetTopMargin(0.15);
        bottomPad->SetTopMargin(0.04);
        bottomPad->SetBottomMargin(0.25);
        bottomPad->SetGridy();
        topPad   ->SetLeftMargin(0.15);
        topPad   ->SetRightMargin(0.1);
        bottomPad->SetLeftMargin(0.15);
        bottomPad->SetRightMargin(0.1);
        topPad   ->Draw();
        bottomPad->Draw();
        topPad   ->cd();

        // Set the Titles ----
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_universe.at(k_up).at(k)->SetTitle(var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_eff)
            h_universe.at(k_up).at(k)->SetTitle(var_labels_eff.at(var).c_str());
        else
            h_universe.at(k_up).at(k)->SetTitle(var_labels_events.at(var).c_str());
        
        h_universe.at(k_up).at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));
        h_universe.at(k_up).at(k)->GetXaxis()->SetTitle("");
        h_universe.at(k_up).at(k)->GetXaxis()->SetLabelSize(0);
        
        
        // Draw the histograms --
        h_universe.at(k_up).at(k)->Draw("hist");
        cv_hist_vec.at(var).at(k)->Draw("hist,same");
        
        if (!single_var)
            h_universe.at(k_dn).at(k)->Draw("hist,same");

        c->Update();

        // Now find the max bin val so we can set the scale of the histogram correctly ----
        double scale_val = h_universe.at(k_up).at(k)->GetMaximum();
        
        if (scale_val < cv_hist_vec.at(var).at(k)->GetMaximum())
            scale_val = cv_hist_vec.at(var).at(k)->GetMaximum();
        
        if (scale_val <  h_universe.at(k_dn).at(k)->GetMaximum())
            scale_val =  h_universe.at(k_dn).at(k)->GetMaximum();

        h_universe.at(k_up).at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        // Fixed scaling for differential cross section
        if (vars.at(var) != "integrated"){
            // h_universe.at(k_up).at(k)->GetYaxis()->SetRangeUser(0, 0.5e-39);
        }

        // Now create the legend, type of legend depends on whether we got a singular variaion or and up/dn type
        TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        
        if (!single_var){
            leg->AddEntry(h_universe.at(k_up).at(k), Form("%s +1 #sigma", label_pretty.c_str()), "l");
            leg->AddEntry(cv_hist_vec.at(var).at(k),         "CV", "l");
            leg->AddEntry(h_universe.at(k_dn).at(k), Form("%s -1 #sigma", label_pretty.c_str()), "l");
        }
        else {
            leg->AddEntry(h_universe.at(k_up).at(k), Form("%s", label_pretty.c_str()), "l");
            leg->AddEntry(cv_hist_vec.at(var).at(k),         "CV", "l");
        }
        
        leg->Draw();

        // Now lets make the ratio plot (actiall percent difference here, but same idea) ----
        bottomPad->cd();
        
        // Up % diff to CV ----
        TH1D* h_err_up = (TH1D *) h_universe.at(k_up).at(k)->Clone("h_ratio_up");
        h_err_up->Add(cv_hist_vec.at(var).at(k), -1);
        h_err_up->Divide(cv_hist_vec.at(var).at(k));
        h_err_up->SetLineWidth(2);
        h_err_up->SetLineColor(kGreen+2);
        h_err_up->Scale(100);

        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_err_up->SetTitle(var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_eff)
            h_err_up->SetTitle(var_labels_eff.at(var).c_str());
        else
            h_err_up->SetTitle(var_labels_events.at(var).c_str());

        
        // Dn % diff to CV ----
        TH1D* h_err_dn = (TH1D *)h_universe.at(k_dn).at(k)->Clone("h_ratio_dn");
        h_err_dn->Add(cv_hist_vec.at(var).at(k), -1);
        h_err_dn->Divide(cv_hist_vec.at(var).at(k));
        h_err_dn->SetLineWidth(2);
        h_err_dn->SetLineColor(kRed+2);
        h_err_dn->Scale(100);

        // Customise a bit of the options
        h_err_up->GetXaxis()->SetLabelSize(0.13);
        h_err_up->GetXaxis()->SetTitleOffset(0.9);
        h_err_up->GetXaxis()->SetTitleSize(0.13);
        h_err_up->GetYaxis()->SetLabelSize(0.13);
        h_err_up->GetYaxis()->SetRangeUser(-50, 50);
        h_err_up->GetYaxis()->SetTitleSize(12);
        h_err_up->GetYaxis()->SetTitleFont(44);
        h_err_up->GetYaxis()->CenterTitle();
        h_err_up->GetYaxis()->SetTitleOffset(1.5);
        h_err_up->SetTitle(" ");
        h_err_up->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        h_err_up->GetYaxis()->SetTitle("\% change from CV");
        h_err_up->Draw("hist,same");

        if (!single_var)
            h_err_dn->Draw("hist,same");

        bottomPad->Update();

        // Draw on the percentage errors manually so they dont overlap
        TLatex* text_up, *text_dn;
        for (int bin = 1; bin < h_err_up->GetNbinsX()+1; bin++){
            
            double bin_up_max = h_err_up->GetBinContent(bin);
            double bin_dn_max = h_err_dn->GetBinContent(bin);

            double shift_up = 0;
            if (bin_up_max >= 0) shift_up = 6;
            else shift_up = -10;

            double shift_dn = 0;
            if (bin_dn_max >= 0) shift_dn = 6;
            else shift_dn = -10;

            if (!single_var){
                if (shift_up < 0 && shift_dn < 0) shift_dn-= 10;
                if (shift_up >= 0 && shift_dn >= 0) shift_up+= 10;
            }
            
            text_up = new TLatex(h_err_up->GetXaxis()->GetBinCenter(bin), bin_up_max+shift_up, Form("%4.1f", h_err_up->GetBinContent(bin)));
            text_dn = new TLatex(h_err_dn->GetXaxis()->GetBinCenter(bin), bin_dn_max+shift_dn, Form("%4.1f", h_err_dn->GetBinContent(bin)));
            text_up->SetTextAlign(21);
            text_up->SetTextColor(kGreen+2);
            text_up->SetTextFont(gStyle->GetTextFont());
            text_up->SetTextSize(0.07);
            text_up->Draw();
            text_dn->SetTextAlign(21);
            text_dn->SetTextColor(kRed+2);
            text_dn->SetTextFont(gStyle->GetTextFont());
            text_dn->SetTextSize(0.07);
            if (!single_var) text_dn->Draw();

        }
        
        // Print the histograms
        c->Print(Form("plots/%s_%s_%s.pdf",sys_label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_err_up;
        delete h_err_dn;
    }

}
// -----------------------------------------------------------------------------
/*
    // Function to make a multisim plot:

    std::string sys_label - The name of the systematic in the file. Choices are: 
    
    weightsGenie, weightsReint, weightsPPFX, MCStats
    ---
    int var - the variable to plot the distributions as a function of. Choices are:
    
    k_integrated, k_reco_el_E, k_true_el_E
    ---
    std::string label_pretty -- name that goes in the legend, choose something appropriate:
    e.g if plotting GENIE, then use GENIE All etc.
    ---
    int universes -- the number of universes. You can use (or lower):
    weightsGenie : 500
    weightsReint : 1000
    weightsPPFX  : 500
    MCStats      : 1000
    ---
    TFile* f -- the input TFile
*/
void PlotMultisim(std::string label, int var, std::string label_pretty, int universes, TFile *f){
    
    TCanvas *c;
    TPad *topPad;
    TPad *bottomPad;

    std::vector<std::vector<TH1D*>> h_universe; // Universe, <gen/sig/xsec etc>
    std::vector<std::vector<TH1D*>> h_err;      // Vector to store the errors

    // Set the x Bins
    std::vector<double> bins;
    if (vars.at(var) == "integrated")
        bins  = { 0.0, 1.1 };
    else
        bins = { 0.0, 0.23, 0.41, 0.65, 0.94, 1.35, 1.87, 2.32, 6.0};
    
    // Get the number of bins and the right vector
    int const nbins = bins.size()-1;
    double* edges = &bins[0]; // Cast to an array 

    // Create the fancier looking 2D histograms
    std::vector<TH2D*> h_universe_2D;
    h_universe_2D.resize(xsec_types.size());
    
    // Set 2D bins for integrated bins
    if (vars.at(var) == "integrated"){
        h_universe_2D.at(k_xsec_sel)           = new TH2D("h_2D_sel",      "", nbins, edges, 100, 2000, 3000);
        h_universe_2D.at(k_xsec_bkg)           = new TH2D("h_2D_bkg",      "", nbins, edges, 100, 200, 1000);
        h_universe_2D.at(k_xsec_gen)           = new TH2D("h_2D_gen",      "", nbins, edges, 100, 4000, 12000);
        h_universe_2D.at(k_xsec_sig)           = new TH2D("h_2D_sig",      "", nbins, edges, 100, 1000, 3000);
        h_universe_2D.at(k_xsec_eff)           = new TH2D("h_2D_eff",      "", nbins, edges, 100, 0.15, 0.3);
        h_universe_2D.at(k_xsec_ext)           = new TH2D("h_2D_ext",      "", nbins, edges, 10, 0, 10);
        h_universe_2D.at(k_xsec_dirt)          = new TH2D("h_2D_dirt",     "", nbins, edges, 15, 0, 15);
        h_universe_2D.at(k_xsec_data)          = new TH2D("h_2D_data",     "", nbins, edges, 80, 0, 140);
        h_universe_2D.at(k_xsec_mcxsec)        = new TH2D("h_2D_mcxsec",   "", nbins, edges, 100, 0.5e-39, 3.0e-39);
        h_universe_2D.at(k_xsec_dataxsec)      = new TH2D("h_2D_dataxsec", "", nbins, edges, 100, 0.5e-39, 3.0e-39);
        
    }
    // Set 2D bins for other
    else {
        h_universe_2D.at(k_xsec_sel)           = new TH2D("h_2D_sel",      "", nbins, edges, 100, 0, 2500);
        h_universe_2D.at(k_xsec_bkg)           = new TH2D("h_2D_bkg",      "", nbins, edges, 100, 0, 1200);
        h_universe_2D.at(k_xsec_gen)           = new TH2D("h_2D_gen",      "", nbins, edges, 100, 0, 16000);
        h_universe_2D.at(k_xsec_sig)           = new TH2D("h_2D_sig",      "", nbins, edges, 100, 0, 2700);
        h_universe_2D.at(k_xsec_eff)           = new TH2D("h_2D_eff",      "", nbins, edges, 100, 0, 0.5);
        h_universe_2D.at(k_xsec_ext)           = new TH2D("h_2D_ext",      "", nbins, edges, 30, 0, 30);
        h_universe_2D.at(k_xsec_dirt)          = new TH2D("h_2D_dirt",     "", nbins, edges, 40, 0, 40);
        h_universe_2D.at(k_xsec_data)          = new TH2D("h_2D_data",     "", nbins, edges, 100, 0, 250);
        h_universe_2D.at(k_xsec_mcxsec)        = new TH2D("h_2D_mcxsec",   "", nbins, edges, 100, 0, 13.0e-39);
        h_universe_2D.at(k_xsec_dataxsec)      = new TH2D("h_2D_dataxsec", "", nbins, edges, 100, 0, 13.0e-39);
    }
    
    
    // Clone the CV histograms so we can work with them without changing them
    std::vector<TH1D*> cv_hist_vec_clone;
    cv_hist_vec_clone.resize(xsec_types.size());
    
    for (unsigned int h_index = 0; h_index < cv_hist_vec_clone.size(); h_index++){
        
        // Clone the histogram
        cv_hist_vec_clone.at(h_index) = (TH1D*)cv_hist_vec.at(var).at(h_index)->Clone(Form("h_%s_clone", xsec_types.at(h_index).c_str() ));
        
        // Customise
        cv_hist_vec_clone.at(h_index)->SetLineWidth(2);
        cv_hist_vec_clone.at(h_index)->SetLineColor(kBlack);
    }

    // Now lets get the multisim histograms from the file

    // Resize to the number of universes
    h_universe.resize(universes);
    
    // Resize each universe
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        h_universe.at(uni).resize(xsec_types.size());
    }
    
    // Get the histograms and customise a bit
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            h_universe.at(uni).at(k) = (TH1D*) f->Get(Form( "%s/%s/h_run1_%s_%i_%s_%s", label.c_str(), vars.at(var).c_str(), label.c_str(), uni ,vars.at(var).c_str(), xsec_types.at(k).c_str()));

            // Customise
            h_universe.at(uni).at(k)->SetLineWidth(1);
            h_universe.at(uni).at(k)->SetLineColor(kAzure+5);

            // Loop over the bins and fill the 2D histogram
            for (int bin = 0; bin < h_universe.at(uni).at(k)->GetNbinsX(); bin++){
                h_universe_2D.at(k)->Fill( bins.at(bin), h_universe.at(uni).at(k)->GetBinContent(bin+1));
            }
            
        }
    }

    // Create vector of systematic error for each histogram
    std::vector<std::vector<double>> sys_err;
    sys_err.resize(xsec_types.size());
    
    // We now want to get the standard deviation of all universes wrt to the cv
    for (unsigned int i = 0; i < sys_err.size(); i++){
        sys_err.at(i).resize(cv_hist_vec.at(var).at(0)->GetNbinsX(), 0.0);
    }

    // Loop over universes
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over histograms
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

            // Loop over the bins
            for (int bin = 1; bin < cv_hist_vec.at(var).at(0)->GetNbinsX()+1; bin++){
                double uni_x_content = h_universe.at(uni).at(k)->GetBinContent(bin);
                double cv_x_content  = cv_hist_vec.at(var).at(k)->GetBinContent(bin);

                sys_err.at(k).at(bin-1) += ( uni_x_content - cv_x_content) * ( uni_x_content - cv_x_content);
            }
            
        }
    }

    // Loop over the histograms
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        
        // loop over the bins
        for (unsigned int bin = 0; bin < sys_err.at(k).size(); bin ++){
            sys_err.at(k).at(bin) = std::sqrt(sys_err.at(k).at(bin) / h_universe.size());
        }
    }

    // Now we have the systematic error computed we should now set the bin error of the CV clone to the std
    // Loop over universes
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over histograms
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

            // Loop over the bins
            for (int bin = 1; bin < cv_hist_vec.at(var).at(0)->GetNbinsX()+1; bin++){
                cv_hist_vec_clone.at(k)->SetBinError(bin, sys_err.at(k).at(bin-1));
            }
            
        }
    }

    // Now we want to draw them
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

        c         = new TCanvas("c", "c", 500, 500);
        topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        
        // Set options for the TPad and draw them
        topPad   ->SetBottomMargin(0.05);
        topPad   ->SetTopMargin(0.15);
        bottomPad->SetTopMargin(0.04);
        bottomPad->SetBottomMargin(0.25);
        bottomPad->SetGridy();
        topPad   ->SetLeftMargin(0.15);
        topPad   ->SetRightMargin(0.1);
        bottomPad->SetLeftMargin(0.15);
        bottomPad->SetRightMargin(0.1);
        topPad   ->Draw();
        bottomPad->Draw();
        topPad   ->cd();

        // Set the title
        h_universe.at(0).at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));

        // Scale the histogram to the right value
        double scale_val = h_universe.at(0).at(k)->GetMaximum();
        
        if (scale_val <  h_universe.at(k_dn).at(k)->GetMaximum())
            scale_val =  h_universe.at(k_dn).at(k)->GetMaximum();

        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_universe_2D.at(k)->SetTitle(var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_eff)
            h_universe_2D.at(k)->SetTitle(var_labels_eff.at(var).c_str());
        else
            h_universe_2D.at(k)->SetTitle(var_labels_events.at(var).c_str());


        // Draw and customise
        h_universe_2D.at(k)->Draw("colz,same");
        h_universe_2D.at(k)->GetXaxis()->SetTitle("");
        h_universe_2D.at(k)->GetXaxis()->SetLabelSize(0);
        h_universe_2D.at(k)->SetTitle(xsec_types_pretty.at(k).c_str());
        h_universe_2D.at(k)->GetYaxis()->SetTitleSize(0.04);
        h_universe_2D.at(k)->GetYaxis()->SetLabelSize(0.05);

        if (xsec_types.at(k) != "data_xsec") {
            cv_hist_vec_clone.at(k)->SetLineColor(kRed+1);
            cv_hist_vec_clone.at(k)->SetLineWidth(2);
            cv_hist_vec_clone.at(k)->SetLineStyle(7);
            cv_hist_vec_clone.at(k)->SetFillStyle(0);
            cv_hist_vec_clone.at(k)->SetMarkerColor(kRed+1);
            cv_hist_vec_clone.at(k)->Draw("E2,hist,same");
        }
        else { 
            cv_hist_vec_clone.at(k)->SetLineColor(kRed+1);
            cv_hist_vec_clone.at(k)->Draw("E,same");
        }

        h_universe.at(0).at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        // Now lets plot the legend ----
        TLegend *leg;
        if (xsec_types.at(k)!= "eff")
            leg = new TLegend(0.5, 0.65, 0.85, 0.8);
        else    
            leg = new TLegend(0.5, 0.1, 0.85, 0.35);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_universe_2D.at(k), Form("%s", label_pretty.c_str()), "");
        
        if (xsec_types.at(k) != "data_xsec")
            leg->AddEntry(cv_hist_vec_clone.at(k),           "CV (Sys Only)", "f");
        else
            leg->AddEntry(cv_hist_vec_clone.at(k),           "CV (Sys Only)", "le");
        
        leg->Draw();


        // Now lets draw the uncertainty on the bottom pad
        bottomPad->cd();

        // error is STD from the CV
        TH1D* h_err = (TH1D *)cv_hist_vec_clone.at(k)->Clone("h_err");
        
        // Loop over the bins in the up error, and set the bin content to be the percent difference
        for (int g = 1; g < h_err->GetNbinsX()+1; g++){
            h_err->SetBinContent(g, 100 * h_err->GetBinError(g)/h_err->GetBinContent(g));
        }

        // Customise the histogram
        h_err->SetLineWidth(2);
        h_err->SetLineColor(kRed+1);
        h_err->GetYaxis()->SetRangeUser(0, 50);
        h_err->GetXaxis()->SetLabelSize(0.13);
        h_err->GetXaxis()->SetTitleOffset(0.9);
        h_err->GetXaxis()->SetTitleSize(0.13);
        h_err->GetYaxis()->SetLabelSize(0.13);
        h_err->GetYaxis()->SetNdivisions(4, 0, 0, kTRUE);
        h_err->GetYaxis()->SetTitleSize(12);
        h_err->GetYaxis()->SetTitleFont(44);
        h_err->GetYaxis()->CenterTitle();
        h_err->GetYaxis()->SetTitleOffset(1.5);

        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_err->SetTitle(var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_eff)
            h_err->SetTitle(var_labels_eff.at(var).c_str());
        else
            h_err->SetTitle(var_labels_events.at(var).c_str());

        h_err->SetTitle(" ");
        h_err->SetMarkerColor(kBlack);
        h_err->SetLineStyle(1);
        h_err->SetLineColor(kBlack);
        h_err->GetYaxis()->SetTitle("\% Uncertainty");
        h_err->SetMarkerSize(4);
        
        if (k != k_xsec_dirt && k != k_xsec_ext)
            h_err->Draw("hist, text00");
        
        gStyle->SetPaintTextFormat("4.1f");
        gStyle->SetPalette(56);
        gStyle->SetPalette(kBlueGreenYellow);

        // Save the histogram
        c->Print(Form("plots/%s_%s_%s.pdf", label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_universe_2D.at(k);

    }

}
// -----------------------------------------------------------------------------
/*
    // Function to make a detvar unisim plot:

    std::string label - The name of the systematic in the file. Choices are: 
    
    LYRayleigh, LYDown, SCE, Recomb2, WireModX, WireModYZ, WireModThetaXZ, WireModThetaYZ_withSigmaSplines
    
    ---
    int var - the variable to plot the distributions as a function of. Choices are:
    
    k_integrated, k_reco_el_E, k_true_el_E
    ---
    int detvar_index -- The enum of the detector variation. Choices are:
    k_LYRayleigh, k_LYDown, k_SCE, k_Recomb2, k_WireModX, k_WireModYZ, k_WireModThetaXZ, k_WireModThetaYZ_withSigmaSplines

    ---
    std::string label_pretty -- name that goes in the legend, we have a vector that can help us with this and use the enum like this:
    
    e.g detvar_string_pretty.at(k_WireModYZ) or detvar_string_pretty.at(k_Recomb2) etc.
    
    ---
    TFile* f -- the input TFile
*/
void PlotDetVar(std::string label, int var, int detvar_index, std::string label_pretty, TFile *f){

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c;

    std::vector<TH1D*> h_universe;
    std::vector<TH1D*> h_CV;
    
    // Resize to the number cross section types e.g. bkg,eff etc.
    h_universe.resize(xsec_types.size());
    h_CV.resize(xsec_types.size());

    TH1D* h_temp;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < h_universe.size(); k++){
        
        // Get the universe histograms
        h_temp = (TH1D*)f->Get(Form( "%s/%s/h_run1_CV_0_%s_%s", label.c_str(), vars.at(var).c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));
        h_universe.at(k) = (TH1D*)h_temp->Clone();

        // Scale the variation to the CV POT
        double scale_fact = POT_v.at(k_CV) / POT_v.at(detvar_index);
        
        // Scale the histograms, but only in the case of MC variables
        if (xsec_types.at(k) != "ext" && xsec_types.at(k) != "dirt" && xsec_types.at(k) != "data")
            h_universe.at(k)->Scale(scale_fact);

        // Get the CV histograms
        h_CV.at(k) = (TH1D*)f->Get(Form( "detvar_CV/%s/h_run1_CV_0_%s_%s", vars.at(var).c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k)->SetLineWidth(2);
        h_universe.at(k)->SetLineColor(kGreen+2);
        h_universe.at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        // Customise
        h_CV.at(k)->SetLineWidth(2);
        h_CV.at(k)->SetLineColor(kBlack);
        h_universe.at(k)->SetLineColor(kGreen+2);
    }

    // Now we want to draw them
    for (unsigned int k = 0; k < h_universe.size(); k++){
        c = new TCanvas("c", "c", 500, 500);
        topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        
        // Set options for the TPad and draw them
        topPad   ->SetBottomMargin(0.05);
        topPad   ->SetTopMargin(0.15);
        bottomPad->SetTopMargin(0.04);
        bottomPad->SetBottomMargin(0.25);
        bottomPad->SetGridy();
        topPad   ->SetLeftMargin(0.15);
        topPad   ->SetRightMargin(0.1);
        bottomPad->SetLeftMargin(0.15);
        bottomPad->SetRightMargin(0.1);
        topPad   ->Draw();
        bottomPad->Draw();
        topPad   ->cd();
        
        // Set the title and draw
        h_universe.at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));
        h_universe.at(k)->GetXaxis()->SetTitle("");
        h_universe.at(k)->GetXaxis()->SetLabelSize(0);

        h_universe.at(k)->Draw("hist");
        h_CV.at(k)->Draw("hist,same");

        c->Update();

        // Scale the histogram
        double scale_val = h_universe.at(k)->GetMaximum();
        if (scale_val < h_CV.at(k)->GetMaximum()) scale_val = h_CV.at(k)->GetMaximum();

        h_universe.at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        // Draw the TLegend
        TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_universe.at(k), label_pretty.c_str(), "l");
        leg->AddEntry(h_CV.at(k), "CV", "l");
        leg->Draw();

        // Now lets do the bottom pad

        bottomPad->cd();
        
        // Up ratio to CV
        TH1D* h_err = (TH1D *)h_universe.at(k)->Clone("h_ratio");
        h_err->Add(h_CV.at(k), -1);
        h_err->Divide(h_CV.at(k));
        h_err->SetLineWidth(2);
        h_err->SetLineColor(kGreen+2);
        h_err->Scale(100);
        
        // Set the Titles
        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec)
            h_err->SetTitle(var_labels_xsec.at(var).c_str());
        else if (k == k_xsec_eff)
            h_err->SetTitle(var_labels_eff.at(var).c_str());
        else
            h_err->SetTitle(var_labels_events.at(var).c_str());
        
        // Set some options for the ratio histograms
        h_err->GetXaxis()->SetLabelSize(0.13);
        h_err->GetXaxis()->SetTitleOffset(0.9);
        h_err->GetXaxis()->SetTitleSize(0.13);
        h_err->GetYaxis()->SetLabelSize(0.13);
        h_err->GetYaxis()->SetRangeUser(-50, 50);
        h_err->GetYaxis()->SetTitleSize(12);
        h_err->GetYaxis()->SetTitleFont(44);
        h_err->GetYaxis()->CenterTitle();
        h_err->GetYaxis()->SetTitleOffset(1.5);
        h_err->SetTitle(" ");
        h_err->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        h_err->GetYaxis()->SetTitle("\% change from CV");
        h_err->Draw("hist,same");

        bottomPad->Update();

        // Draw on the percentage errors manually so they dont overlap
        TLatex* text_up, *text_dn;
        for (int bin = 1; bin < h_err->GetNbinsX()+1; bin++){
            
            double bin_up_max = h_err->GetBinContent(bin);

            double shift_up = 0;
            if (bin_up_max >= 0) shift_up = 6;
            else shift_up = -10;
            
            text_up = new TLatex(h_err->GetXaxis()->GetBinCenter(bin), bin_up_max+shift_up, Form("%4.1f", h_err->GetBinContent(bin)));
            text_up->SetTextAlign(21);
            text_up->SetTextColor(kGreen+2);
            text_up->SetTextFont(gStyle->GetTextFont());
            text_up->SetTextSize(0.07);
            text_up->Draw();
        }
        
        // Save the histgoram to file
        c->Print(Form("plots/%s_%s_%s.pdf", label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_err;
    }
}