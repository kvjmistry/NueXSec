/* 

Macro to Make the weight histograms for the Central value weighting

This can and should be extended for making the systematic variations

USAGE: root -l -b -q 'Weight_Histograms.C("fhc/rhc")'

*/
// ------------------------------------------------------------------------------------------------------------
bool GetHist(TFile* f, TH1D* &h, TString string){
    h = (TH1D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool GetHist(TFile* f, TH2D* &h, TString string){
    h = (TH2D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool GetFile(TFile* &f , TString string){
    f = TFile::Open(string);
    
    if (f == NULL) {
        std::cout << "failed to get:\t" << string << "\tThis file might not exist in the file" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
//  Function to calculate histogram ratios and return them as a vector of Th2D
void CalcRatio(TFile *fin, TFile *fout,int flav, const char* mode){
    TH2D *h_CV_UW, * h_CV_PPFX;
    bool boolhist;
    
    fin->cd();
    

    
    std::cout << "Calculating Histograms...\t" << mode << std::endl;

    // 2D
    //  Get the histograms
    boolhist = GetHist(fin, h_CV_UW, Form("%s/Detsmear/%s_unweighted_AV_TPC_2D",   mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the unweighted CV
    boolhist = GetHist(fin, h_CV_PPFX, Form("%s/Detsmear/%s_CV_AV_TPC_2D", mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the PPFX weighted CV

    TH2D *hratio_2D = (TH2D*) h_CV_PPFX->Clone(Form("h_2D_CV_UW_PPFX_ratio_%s",mode));

    hratio_2D->Divide(h_CV_UW); // Divide hists

    hratio_2D->SetOption("colz");

    fout->cd();
    hratio_2D->Write("",TObject::kOverwrite);

    // 1D
    fin->cd();
    TH1D *h_1D_CV_UW, * h_1D_CV_PPFX;

    boolhist = GetHist(fin, h_1D_CV_UW, Form("%s/Detsmear/%s_UW_AV_TPC",   mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the unweighted CV
    boolhist = GetHist(fin, h_1D_CV_PPFX, Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the PPFX weighted CV

    TH1D *hratio = (TH1D*) h_1D_CV_PPFX->Clone(Form("h_1D_CV_UW_PPFX_ratio_%s",mode));

    hratio->Divide(h_1D_CV_UW); // Divide hists

    hratio->SetOption("E");

    fout->cd();
    hratio->Write("",TObject::kOverwrite);

    // Theta

    fin->cd();
    TH1D *h_1D_CV_UW_theta, * h_1D_CV_PPFX_theta;

    boolhist = GetHist(fin, h_1D_CV_UW_theta, Form("%s/Detsmear/Th_%s_UW_TPC",   mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the unweighted CV
    boolhist = GetHist(fin, h_1D_CV_PPFX_theta, Form("%s/Detsmear/Th_%s_CV_TPC", mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the PPFX weighted CV

    TH1D *hratio_theta = (TH1D*) h_1D_CV_PPFX_theta->Clone(Form("h_1D_theta_CV_UW_PPFX_ratio_%s",mode));

    hratio_theta->Divide(h_1D_CV_UW_theta); // Divide hists

    hratio_theta->SetOption("E");

    fout->cd();
    hratio_theta->Write("",TObject::kOverwrite);

}
// -----------------------------------------------------------------------------
void Weight_Histograms(std::string horn){

    enum flavours {
        knue,
        knuebar,
        knumu,
        knumubar,
        k_flav_MAX
    };

    std::vector<const char*> modes = {"nue", "nuebar", "numu", "numubar"};


    std::vector<TH2D*> hists;
    hists.resize(k_flav_MAX);

    // File with CV
    TFile *fin;

    bool boolfile;
    
    if (horn == "fhc" )boolfile = GetFile(fin , "output_fhc_uboone_run0_merged.root");
    else if (horn == "rhc") boolfile = GetFile(fin , "output_rhc_uboone_run0_merged.root"); 
    else {
        std::cout << "You misepelt the horn config...."<< std::endl;
        gSystem->Exit(0);
    }
    
    if (boolfile == false) gSystem->Exit(0); // Most up to date version of CV

    TFile *fout;
    
    if      (horn == "fhc") fout = new TFile("f_flux_CV_weights_fhc.root", "RECREATE");
    else if (horn == "rhc") fout = new TFile("f_flux_CV_weights_rhc.root", "RECREATE");

    for (unsigned int k =0 ; k < k_flav_MAX; k++){
        CalcRatio(fin, fout,  k, modes.at(k));
    }

}
// -----------------------------------------------------------------------------
