// Script to optimise the bin widths in the truth electron distribution

void GetFitResult(double &mean, double &sigma, float bin_lower_edge, float bin_upper_edge, TTree* tree, bool save_hist, bool &converged){
    
    TCut generic_query = "(classifcation.c_str()==\"nue_cc\" || classifcation.c_str()==\"nuebar_cc\") && !gen"; // This gets selected signal events
    TCut bin_query = Form("shr_energy_tot_cali > %f && shr_energy_tot_cali < %f", bin_lower_edge, bin_upper_edge);
    
    TCanvas *c = new TCanvas();

    // Draw the Query
    tree->Draw("elec_e", generic_query && bin_query);
    
    // Get the histogram from the pad
    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
    
    // Fit it with a Gaussian
    htemp->Fit("gaus");
    
    // Get the fit result
    TF1 *fit_gaus = htemp->GetFunction("gaus");
    
    // Draw the histogram
    if (save_hist) {
        htemp->SetLineWidth(2);
        htemp->SetLineColor(kBlack);
    }
    htemp->Draw("hist");
    fit_gaus->Draw("same");

    mean  = fit_gaus->GetParameter(1);
	sigma = fit_gaus->GetParameter(2);

    if (sigma*2 >= bin_upper_edge - bin_lower_edge - 0.01 && sigma*2 <= bin_upper_edge - bin_lower_edge + 0.01){
        std::cout << "Fit has converged!: " << 2*sigma/(bin_upper_edge - bin_lower_edge) << std::endl;
        converged = true;
    }

    if (save_hist) c->Print(Form("plots/bins_%0.2fGeV_to_%0.2f_GeV.pdf",bin_lower_edge, bin_upper_edge ));


}


void optmise_bins(){


    // Load in the tfile and tree
    TFile *f = TFile::Open("/Users/kvjmistry/Documents/work/NueXSec/Analysis/files/trees/nuexsec_selected_tree_mc_run1.root");
    TTree *tree = (TTree*) f->Get("mc_tree");


    double mean{0.0}, sigma{0.0};
    bool converged = false;

    // call function which draws the tree to a canvas, fits the tree and returns the fit parameter
    // If the fit has 2xSTD = the reco bin size then we have successfully optimised the bin
    float lower_bin = 0.25;
    
    // Loop over the bins
    for (float bin = 0; bin < 4; bin++ ){
        std::cout << "\n\033[0;34mTrying to optimise the next bin\033[0m\n"<< std::endl;
        converged = false;

        // Slide upper bin value till we get 2xthe STD of the fit
        for (float i = lower_bin+0.1; i < 4.0; i+=0.01) {
            std::cout << "\n\033[0;34mTrying Bin: " << i << "GeV\033[0m\n"<< std::endl;
            GetFitResult(mean, sigma, lower_bin, i, tree, false, converged);

            // If it converged, do it again and print the canvas then break
            if (converged) {
                GetFitResult(mean, sigma, lower_bin, i, tree, true, converged);
                std::cout << "\n\033[0;34mMean: " << mean << "  Sigma: " << sigma<< "\033[0m\n"<< std::endl;
                
                // Reset the lower bin value
                lower_bin = i;
                break;
            }

        }
    }
    
}


