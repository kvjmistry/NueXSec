// -----------------------------------------------------------------------------
double GetNuMIAngle(double px, double py, double pz, std::string direction){

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
    if (direction == "beam") beamdir = {0 , 0 , 1};
    
    // Get the angle wrt to the target to detector direction
    else if (direction == "target") {
        beamdir = {5502, 7259, 67270};
        beamdir = beamdir.Unit(); // Get the direction
    }
    else {
        std::cout << "Warning unknown angle type specified, you should check this" << std::endl;
    }
    
    double angle = BeamCoords.Angle(beamdir) * 180 / 3.1415926;


    // Create vectors to get the angle in the yz and xz planes
    TVector3 BeamCoords_yz = { 0, BeamCoords.Y(), BeamCoords.Z() }; // Angle upwards
    TVector3 BeamCoords_xz = { BeamCoords.X(), 0, BeamCoords.Z() }; // Angle across

    // if (theta > 50 ) std::cout <<"Theta: " << theta << "   UP: " << BeamCoords_yz.Angle(beam_dir) * 180 / 3.1415926 << "  Across: " << BeamCoords_xz.Angle(beam_dir) * 180 / 3.1415926 << std::endl;

    // std::cout << theta << std::endl;

    return angle;
}



void plot_ang_baseline(){


    TFile *f = TFile::Open("~/Desktop/neutrinoselection_filt_run1_overlay_intrinsic.root", "READ");

    TTree *t = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");

    double baseline, par_decay_vz;
    float true_nu_px, true_nu_py, true_nu_pz;
    t->SetBranchAddress("baseline", &baseline);
    t->SetBranchAddress("par_decay_vz", &par_decay_vz);
    t->SetBranchAddress("true_nu_px", &true_nu_px);
    t->SetBranchAddress("true_nu_py", &true_nu_py);
    t->SetBranchAddress("true_nu_pz", &true_nu_pz);

    int  tree_total_entries = t->GetEntries();

    TH2D* h_baseline = new TH2D("h_baseline", ";Baseline [m]; Angle [deg]", 100, 50, 710, 100, 0, 180);
    TH2D* h_zpos = new TH2D("h_zpos", ";z pos Beamline [m]; Angle [deg]", 100, 0, 1000, 180, 0, 180);

    // Event loop
    for (int ievent = 0; ievent < tree_total_entries; ievent++){
        
        // Alert the user
        if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;

        // Get the entry in the tree
        t->GetEntry(ievent); 


        double ang = GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "beam");
        h_baseline->Fill(baseline, ang);
        h_zpos->Fill(par_decay_vz/100, ang);

    }

    TCanvas *c = new TCanvas("", "", 500, 500);
    h_baseline->Draw("colz");

    TCanvas *c2 = new TCanvas("", "", 500, 500);
    h_zpos->Draw("colz");

    std::cout << 100* h_zpos->Integral(0, h_zpos->GetNbinsX()+1, 7, 12)/ h_zpos->Integral() << std::endl;





}