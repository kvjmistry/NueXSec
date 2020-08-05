#include "../include/Utility.h"

// -----------------------------------------------------------------------------
void Utility::Initalise(const char* variation, bool overwritePOT, const char* run_period){

    std::cout << "Initialising Utility Class..." << std::endl;

    std::string variation_str = variation;

    // Check the variaition mode setting and see if we want to overwrite the POT
    if (overwritePOT){

        // Loop over the POT config names and overwrite the name of the CV MC POT
        for (unsigned int p = 0; p < confignames.size(); p++){
            
            std::string match_name = Form("Run%s_MC_POT", run_period);
            // std::string match_name = "Run1_MC_POT";
            
            // If matched then overwrite the POT config for the MC to the variation
            if (confignames.at(p) == match_name){
                confignames[p] = match_name + "_" + variation_str;
                std::cout << "New MC POT config to search for is: " << confignames.at(p) << std::endl;
            }

        }
    }

    // Get the congigureation parameters
    std::string line;

    config_v.resize(k_config_MAX, 1.0);

    std::string varname;
    std::string value;
    
    // Loop over the config ist
    for (unsigned int p = 0; p < confignames.size(); p++){

        std::ifstream myfile ("config.txt");

        if (myfile.is_open()) {
            
            // std::cout << confignames.at(p) <<  std::endl;

            // Loop over lines in file
            while ( getline (myfile,line) ) {

                std::istringstream ss(line);
                ss >> varname >> value;

                // Found the correct variation file 
                if (varname == confignames.at(p)) {
                    // std::cout << "Found match for: " << varname << " "<< std::stod(value) <<  std::endl;
                    config_v.at(p)= std::stod(value);
                    break;
                }
                
            }
           
       
        }
        else std::cout << "Unable to open file, bad things are going to happen..." << std::endl; 

        myfile.close();
    }

}
// -----------------------------------------------------------------------------
bool Utility::GetFile(TFile* &f, TString string){
    f = TFile::Open(string);
    
    if (f == NULL) {
        std::cout << "failed to get:\t" << string << "\tThis mode wont be used in the selection" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool Utility::GetTree(TFile* f, TTree* &T, TString string){
    T = (TTree*) f->Get(string);
    if (T == NULL) {
        std::cout << "\nfailed to get:\t" << string << "\tThis tree might not exist in the file\n" << std::endl;
        // exit(1);
        return true;
    }
    else {
        return false;
    }
}
// -----------------------------------------------------------------------------
bool Utility::GetHist(TFile* f, TH1D* &h, TString string){
    h = (TH1D*) f->Get(string);
    if (h == NULL) {
        if (verbose) std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool Utility::GetHist(TFile* f, TH2D* &h, TString string){
    h = (TH2D*) f->Get(string);
    if (h == NULL) {
        if (verbose) std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool Utility::GetDirectory(TFile* f, TDirectory* &d, TString string){
    d = (TDirectory*)f->Get(string);
    if (d == NULL) {
        // std::cout << "\nfailed to get:\t" << string << "\tThis directory might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
void Utility::Tabulate(bool inFV, std::string interaction, std::string classification, int type, std::vector<double> &counter_v, double weight) {

    if (type == k_mc){
        
        // Require in FV condition
        if (inFV){
            if (interaction == "nue_cc_qe")  counter_v.at(k_count_nue_cc_qe)  += weight; 
            if (interaction == "nue_cc_res") counter_v.at(k_count_nue_cc_res) += weight;
            if (interaction == "nue_cc_coh") counter_v.at(k_count_nue_cc_dis) += weight;
            if (interaction == "nue_cc_dis") counter_v.at(k_count_nue_cc_coh) += weight;
            if (interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_mec) += weight;

            if (interaction == "nue_bar_cc_qe")  counter_v.at(k_count_nuebar_cc_qe)  += weight; 
            if (interaction == "nue_bar_cc_res") counter_v.at(k_count_nuebar_cc_res) += weight;
            if (interaction == "nue_bar_cc_coh") counter_v.at(k_count_nuebar_cc_dis) += weight;
            if (interaction == "nue_bar_cc_dis") counter_v.at(k_count_nuebar_cc_coh) += weight;
            if (interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_mec) += weight;

            if (interaction == "numu_cc_qe"  )  counter_v.at(k_count_numu_cc_qe)  += weight;
            if (interaction == "numu_cc_res" )  counter_v.at(k_count_numu_cc_res) += weight;
            if (interaction == "numu_cc_coh" )  counter_v.at(k_count_numu_cc_dis) += weight;
            if (interaction == "numu_cc_dis" )  counter_v.at(k_count_numu_cc_coh) += weight;
            if (interaction == "numu_cc_mec" )  counter_v.at(k_count_numu_cc_mec) += weight;

            if (interaction == "numu_bar_cc_qe")  counter_v.at(k_count_numubar_cc_qe)  += weight;
            if (interaction == "numu_bar_cc_res") counter_v.at(k_count_numubar_cc_res) += weight;
            if (interaction == "numu_bar_cc_coh") counter_v.at(k_count_numubar_cc_dis) += weight;
            if (interaction == "numu_bar_cc_dis") counter_v.at(k_count_numubar_cc_coh) += weight;
            if (interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_mec) += weight;

            // These are all the nus, but now in the fv
            if (interaction == "nue_cc_qe" || interaction == "nue_cc_res" || interaction == "nue_cc_coh" || interaction == "nue_cc_dis" || interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_infv) += weight;
            if (interaction == "nue_bar_cc_qe" || interaction == "nue_bar_cc_res" || interaction == "nue_bar_cc_coh" || interaction == "nue_bar_cc_dis" || interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_infv) += weight;
            
            if (interaction == "numu_cc_qe" || interaction == "numu_cc_res" || interaction == "numu_cc_coh" || interaction == "numu_cc_dis" || interaction == "numu_cc_mec") counter_v.at(k_count_numu_cc_infv) += weight;
            if (interaction == "numu_bar_cc_qe" || interaction == "numu_bar_cc_res" || interaction == "numu_bar_cc_coh" || interaction == "numu_bar_cc_dis" || interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_infv) += weight;
        }

        // These are all the nus, but now in the cryostat volume
        if (interaction == "nue_cc_qe" || interaction == "nue_cc_res" || interaction == "nue_cc_coh" || interaction == "nue_cc_dis" || interaction == "nue_cc_mec") counter_v.at(k_count_nue_cc_incryo) += weight;
        if (interaction == "nue_bar_cc_qe" || interaction == "nue_bar_cc_res" || interaction == "nue_bar_cc_coh" || interaction == "nue_bar_cc_dis" || interaction == "nue_bar_cc_mec") counter_v.at(k_count_nuebar_cc_incryo) += weight;
        
        if (interaction == "numu_cc_qe" || interaction == "numu_cc_res" || interaction == "numu_cc_coh" || interaction == "numu_cc_dis" || interaction == "numu_cc_mec") counter_v.at(k_count_numu_cc_incryo) += weight;
        if (interaction == "numu_bar_cc_qe" || interaction == "numu_bar_cc_res" || interaction == "numu_bar_cc_coh" || interaction == "numu_bar_cc_dis" || interaction == "numu_bar_cc_mec") counter_v.at(k_count_numubar_cc_incryo) += weight;
        
        // Classification
        if (classification == "nue_cc")       counter_v.at(k_count_nue_cc)       += weight;
        if (classification == "nuebar_cc")    counter_v.at(k_count_nuebar_cc)    += weight;
        if (classification == "nu_out_fv")    counter_v.at(k_count_nu_out_fv)    += weight;
        if (classification == "nc")           counter_v.at(k_count_nc)           += weight;
        if (classification == "nc_pi0")       counter_v.at(k_count_nc_pi0)       += weight;
        if (classification == "numu_cc")      counter_v.at(k_count_numu_cc)      += weight;
        if (classification == "numu_cc_pi0")  counter_v.at(k_count_numu_cc_pi0)  += weight;
        if (classification == "cosmic")       counter_v.at(k_count_cosmic)       += weight;
        if (classification == "unmatched")    counter_v.at(k_count_unmatched)    += weight;

        // Total selected MC events
        counter_v.at(k_count_total_mc) += weight;
    
    }
    else if (type == k_data) {
        counter_v.at(k_count_data) += weight; // The weight **should** always be 1 for these
    }
    else if (type == k_ext){
        counter_v.at(k_count_ext)  += weight; // The weight **should** always be 1 for these -- maybe not this if we decide to use the BNB stream
    }
    else if (type == k_dirt){
        counter_v.at(k_count_dirt) += weight;
    }
    else {

        std::cout << "unkown type specified!!!  " << __PRETTY_FUNCTION__ << std::endl;
    }
    
}
// -----------------------------------------------------------------------------
double Utility::GetNuMIAngle(double px, double py, double pz, std::string direction){

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
    // NuMI acos(Pz/P)
    else if (direction == "numi_theta"){
        double p = std::sqrt(BeamCoords.X()*BeamCoords.X() +BeamCoords.Y()*BeamCoords.Y() +  BeamCoords.Z()*BeamCoords.Z() );
        return acos(BeamCoords.Z()/p) * (180 / 3.1415);
    }
    else if (direction == "numi_phi"){
        return atan2(BeamCoords.Y(), BeamCoords.X()) * 180 / 3.1415;
    } 
    else {
        std::cout << "Warning unknown angle type specified, you should check this" << std::endl;
    }
    
    double theta = BeamCoords.Angle(beamdir) * 180 / 3.1415926;


    // Create vectors to get the angle in the yz and xz planes
    TVector3 BeamCoords_yz = { 0, BeamCoords.Y(), BeamCoords.Z() }; // Angle upwards
    TVector3 BeamCoords_xz = { BeamCoords.X(), 0, BeamCoords.Z() }; // Angle across

    // if (theta > 50 ) std::cout <<"Theta: " << theta << "   UP: " << BeamCoords_yz.Angle(beam_dir) * 180 / 3.1415926 << "  Across: " << BeamCoords_xz.Angle(beam_dir) * 180 / 3.1415926 << std::endl;

    // std::cout << theta << std::endl;

    return theta;
}
// -----------------------------------------------------------------------------
bool Utility::in_fv(double x, double y, double z){
    
    // Require the vertex is in the boundary
    if ( x   >= config_v.at(k_config_x1) && x <= config_v.at(k_config_x2)
      && y   >= config_v.at(k_config_y1) && y <= config_v.at(k_config_y2)
      && z   >= config_v.at(k_config_z1) && z <= config_v.at(k_config_z2)  ){
        return true;   // pass
    }  
    else return false; // fail
}
// -----------------------------------------------------------------------------