float GetPPFXCVWeight(float nu_e, int nu_pdg, double angle);
double GetNuMIAngle(double px, double py, double pz);


enum flav {
    k_numu,
    k_numubar,
    k_nue,
    k_nuebar,
    k_FLAV_MAX
};

std::vector<std::string> flav_str = {
    "numu",
    "numubar",
    "nue",
    "nuebar"
};

std::vector<TH2D*> hist_ratio(k_FLAV_MAX);
std::vector<TH2D*> hist_ratio_uw(k_FLAV_MAX);

void add_ppfx_cv(){

    TFile * f_flux = TFile::Open("../../../Analysis/Systematics/output_fhc_uboone_run0.root", "READ");

    // Get the flux histograms
    for (unsigned int f = 0; f < flav_str.size(); f++){
        hist_ratio.at(f)      = (TH2D*) f_flux->Get(Form("%s/Detsmear/%s_CV_AV_TPC_2D", flav_str.at(f).c_str(), flav_str.at(f).c_str()));
        hist_ratio_uw.at(f)   = (TH2D*) f_flux->Get(Form("%s/Detsmear/%s_unweighted_AV_TPC_2D", flav_str.at(f).c_str(), flav_str.at(f).c_str()));

        hist_ratio.at(f)->SetDirectory(0);
        hist_ratio_uw.at(f)->SetDirectory(0);

        // Take the ratio
        hist_ratio.at(f)->Divide(hist_ratio_uw.at(f));

    }

    f_flux->Close();

    // Now load in the TTree
    TFile *f = new TFile("neutrinoselection_filt_run1_overlay_nuwro.root","UPDATE");
    TTree *T = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter"); 
    
    float true_nu_px,true_nu_py, true_nu_pz;
    float nu_e;
    int nu_pdg;
    float _ppfx_cv;
    
    TBranch *bpt = T->Branch("_ppfx_cv",&_ppfx_cv,"_ppfx_cv/F");
    T->SetBranchAddress("true_nu_px",&true_nu_px);
    T->SetBranchAddress("true_nu_py",&true_nu_py);
    T->SetBranchAddress("true_nu_pz",&true_nu_pz);
    T->SetBranchAddress("nu_e",&nu_e);
    T->SetBranchAddress("nu_pdg",&nu_pdg);
    

    Long64_t nentries = T->GetEntries();
    
    for (Long64_t i=0;i<nentries;i++){

        // if (i == 10)
        //     break;

        if (i % 100000 == 0) std::cout << "On entry " << i/100000.0 <<"00k " << std::endl;

        T->GetEntry(i);
        
        double angle = GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz);
        
        _ppfx_cv = GetPPFXCVWeight(nu_e, nu_pdg, angle);

        // std::cout << _ppfx_cv << std::endl;

        bpt->Fill();
    }
    
    // T->Print();
    f->cd("nuselection");
    T->Write("NeutrinoSelectionFilter", TObject::kOverwrite);
    
    // delete f;


}


// -----------------------------------------------------------------------------
float GetPPFXCVWeight(float nu_e, int nu_pdg, double nu_theta){
    
    float weight = 1.0;

    double xbin{1.0},ybin{1.0};

    if (nu_pdg == 14) {
        xbin =  hist_ratio.at(k_numu)->GetXaxis()->FindBin(nu_e);
        ybin =  hist_ratio.at(k_numu)->GetYaxis()->FindBin(nu_theta);
        weight =  hist_ratio.at(k_numu)->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == -14) {
        xbin =  hist_ratio.at(k_numubar)->GetXaxis()->FindBin(nu_e);
        ybin = hist_ratio.at(k_numubar)->GetYaxis()->FindBin(nu_theta);
        weight = hist_ratio.at(k_numubar)->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == 12) {
        xbin = hist_ratio.at(k_nue)->GetXaxis()->FindBin(nu_e);
        ybin = hist_ratio.at(k_nue)->GetYaxis()->FindBin(nu_theta);
        weight = hist_ratio.at(k_nue)->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == -12) {
        xbin = hist_ratio.at(k_nuebar)->GetXaxis()->FindBin(nu_e);
        ybin = hist_ratio.at(k_nuebar)->GetYaxis()->FindBin(nu_theta);
        weight = hist_ratio.at(k_nuebar)->GetBinContent(xbin, ybin);
    }

    // Add some catches to remove unphysical weights
    if (std::isinf(weight))      weight = 1.0; 
    if (std::isnan(weight) == 1) weight = 1.0;
    if (weight > 100)            weight = 1.0;

    return weight;
}


double GetNuMIAngle(double px, double py, double pz){

    // Variables
    TRotation RotDet2Beam;             // Rotations
    TVector3  detxyz, BeamCoords;      // Translations
    std::vector<double> rotmatrix;     // Inputs

    // input detector coordinates to translate
    detxyz = {px, py, pz};     

    // From beam to detector rotation matrix
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam
    // RotDet2Beam.Invert(); // Invert back to the beam to det

    // Rotate to beam coords
    BeamCoords = RotDet2Beam * detxyz;

    TVector3 beamdir = {0 , 0 , 1};;
    
    // Get the angle wrt to the beam
    beamdir = {0 , 0 , 1};

    double angle = BeamCoords.Angle(beamdir) * 180 / 3.1415926;


    // Create vectors to get the angle in the yz and xz planes
    TVector3 BeamCoords_yz = { 0, BeamCoords.Y(), BeamCoords.Z() }; // Angle upwards
    TVector3 BeamCoords_xz = { BeamCoords.X(), 0, BeamCoords.Z() }; // Angle across

    // if (theta > 50 ) std::cout <<"Theta: " << theta << "   UP: " << BeamCoords_yz.Angle(beam_dir) * 180 / 3.1415926 << "  Across: " << BeamCoords_xz.Angle(beam_dir) * 180 / 3.1415926 << std::endl;

    // std::cout << theta << std::endl;

    return angle;
}