#include "../include/SystematicsHelper.h"

// -----------------------------------------------------------------------------
void SystematicsHelper::Initialise(Utility _utility){

    std::cout << "Initalising Systematics Helper..." << std::endl;
    _util = _utility;

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

    // Off beam mode to compare bnb and numi off beam samples
    if (std::string(_util.sysmode) == "ext"){
        var_string = { "NuMI", "BNB" };
    }

    // If we choose this mode then we actually want to use a different initialiser
    if (std::string(_util.sysmode) == "reweight"){
        InitialiseReweightingMode();
        return;
    }

    // Get the POT of the variations from the file
    GetPOT();

    // check if current POT values match with the ones in the files
    //if (strcmp(_util.check_pot) == 0){
    //	_util.check_pot = true;
    //}   
 
    // Resize the file vector
    f_vars.resize(k_vars_MAX);

    // Get the variation files
    for (unsigned int l =0; l < var_string.size(); l++){
        
        // Standard variation mode
        if (std::string(_util.sysmode) == "default")  {
            f_vars.at(l) = new TFile( Form("files/nuexsec_run%s_%s_merged.root", _util.run_period, var_string.at(l).c_str() ), "READ");
        }
        // Off beam mode
        else if (std::string(_util.sysmode) == "ext") {
            f_vars.at(l) = new TFile( Form("files/nuexsec_ext_run%s_%s.root", _util.run_period, var_string.at(l).c_str() ), "READ");
        }
        else {
            std::cout << "Error I dont know what mode you have configured..." << std::string(_util.sysmode) << std::endl;
            exit(1);
        }
        
        
        // Maybe add something here to pickup non processed variation
    }

    // Now loop over events and caluclate the cross section
    MakeHistograms();

}

// -----------------------------------------------------------------------------
void SystematicsHelper::GetPOT(){

    std::cout << "Getting the POT for the variation files" << std::endl;

    POT_v.resize(var_string.size(), 1.0);

    std::string line;

    std::string varname;
    std::string value;

    std::string pot_mode = "_MC_POT_"; // default mode

    if (mode == "ext") pot_mode =  "_EXT_trig_";

    std::string POT_run_config = "Run" + std::string(_util.run_period) + pot_mode;
    
    // Loop over the config ist
    for (unsigned int p = 0; p < var_string.size(); p++){

        std::ifstream myfile ("config.txt");

        if (myfile.is_open()) {
            
            // std::cout << var_string.at(p) <<  std::endl;

            // Loop over lines in file
            while ( getline (myfile,line) ) {

                std::istringstream ss(line);
                ss >> varname >> value;

                // Found the correct variation file 
                std::string POT_name_match = POT_run_config + var_string.at(p);

                if (varname == POT_name_match ) {
                    std::cout << "Found match for: " << varname << " "<< std::stod(value) <<  std::endl;
                    POT_v.at(p)= std::stod(value);
                    break;
                }
                
            }
           
        }
        else std::cout << "Unable to open config file, bad things are going to happen..." << std::endl; 

        myfile.close();
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetVariationProperties(TH1D* h, int index){

    if (mode == "default"){
        if (index == k_CV){
            h->SetLineColor(kBlack);
        }
        else if (index == k_bnb_diffusion){
            h->SetLineColor(kAzure-5);
        }
    }
    else {
        
        if (index == k_NuMI){
            h->SetLineColor(kBlack);
        }
        else if (index == k_BNB){
            h->SetLineColor(kAzure-5);
        }

    }
    
    h->SetLineWidth(2);
    h->GetYaxis()->SetTitle("Entries");

}
// -----------------------------------------------------------------------------
void SystematicsHelper::Draw_Area_Norm(TCanvas* c){
    c->cd();

    TPaveText *pt;

    pt = new TPaveText(0.4, 0.916, 0.4, 0.916,"NDC");
    pt->AddText("Area Normalised");
    pt->SetTextColor(kGreen+2);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
    pt->Draw();
}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetRatioOptions(TH1D* hist ){

    hist->GetXaxis()->SetLabelSize(0.13);
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetXaxis()->SetTitleSize(0.13);
    hist->GetYaxis()->SetLabelSize(0.13);
    hist->GetYaxis()->SetRangeUser(-50, 50);
    hist->GetYaxis()->SetTitleSize(12);
    hist->GetYaxis()->SetTitleFont(44);
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleOffset(1.5);
    hist->SetTitle(" ");

}
// -----------------------------------------------------------------------------
void SystematicsHelper::MakeHistograms(){

    // ---- making histograms for the uncertainties
    std::cout << "Making histograms for the detector systematics uncertainties:" << std::endl;

    // Create the root file to save the variation plots
    // this file will be used later to calculate the sys uncertainties
    TFile *file_sys_var = new TFile(Form("plots/run%s/systvar/run%s_sys_var.root", _util.run_period, _util.run_period),"RECREATE");


	/*

    // name of the histogram in the ntuple file
    // write it without the initial "h_"
    // this name is going to be used to name the pdf file
    std::vector<std::string> vec_hist_name = {"reco_vtx_x_sce",
                                              "reco_vtx_y_sce",
                                              "reco_vtx_z_sce",
                                              "reco_flash_time",
                                              "reco_leading_shower_phi",
                                              "reco_leading_shower_theta",
                                              "reco_shower_multiplicity",
                                              "reco_track_multiplicity",
                                              "reco_topological_score",
                                              "reco_shr_tkfit_dedx_y",
                                              "reco_nu_e",
                                              "reco_shower_energy_tot_cali",
                                              "reco_softwaretrig",
                                              "reco_nslice",
                                              "reco_shower_score",
                                              "reco_shr_tkfit_dedx_max",
                                              "reco_shr_tkfit_dedx_max_with_tracks",
                                              "reco_shr_tkfit_dedx_max_no_tracks",
                                              "reco_shower_to_vtx_dist",
                                              "reco_hits_ratio",
                                              "reco_CosmicIPAll3D",
                                              "reco_contained_fraction",
                                              "reco_shrmoliereavg",
                                              "reco_shower_energy_cali",
                                              "reco_flash_pe"};

    // x axis label
    std::vector<std::string> vec_axis_label = {"Reco Vertex X [cm]",
                                               "Reco Vertex Y [cm]",
                                               "Reco Vertex Z [cm]",
                                               "Flash Time [#mus]",
                                               "Leading Shower Phi [degrees]",
                                               "Leading Shower Theta [degrees]",
                                               "Shower Multiplicty",
                                               "Track Multiplicity",
                                               "Topological Score",
                                               "Collection Plane dEdx (track fitter) [MeV/cm]",
                                               "Reconstructed Neutrino Energy [GeV]",
                                               "Reconstructed Leading Shower Energy [GeV]",
                                               "Software Trigger",
                                               "nslice",
                                               "Shower Score",
                                               "reco_shr_tkfit_dedx_max",
                                               "reco_shr_tkfit_dedx_max_with_tracks",
                                               "reco_shr_tkfit_dedx_max_no_tracks",
                                               "reco_shower_to_vtx_dist",
                                               "reco_hits_ratio",
                                               "reco_CosmicIPAll3D",
                                               "reco_contained_fraction",
                                               "reco_shrmoliereavg",
                                               "reco_shower_energy_cali",
                                               "reco_flash_pe"};     */


    for (unsigned int i = 0 ; i < _util.k_cuts_MAX; i++){

        if (mode == "default"){ // default detector systematics mode

            // Create the directory for sysvar
            _util.CreateDirectory("/systvar/comparisons/cuts/" + _util.cut_dirs.at(i));

            // create the directory for detvar
            _util.CreateDirectory("/detvar/comparisons/cuts/" + _util.cut_dirs.at(i));

            for(unsigned int j=0; j < _util.vec_hist_name.size(); j++){

                SysVariations(Form("h_%s", _util.vec_hist_name.at(j).c_str()), Form("plots/run%s/systvar/comparisons/cuts/%s/%s.pdf", _util.run_period, _util.cut_dirs.at(i).c_str(), _util.vec_hist_name.at(j).c_str()),
                            _util.cut_dirs.at(i), _util.vec_axis_label.at(j).c_str(), _util.cut_dirs.at(i).c_str(), _util.vec_hist_name.at(j).c_str(), file_sys_var );

                PlotVariations(Form("h_%s", _util.vec_hist_name.at(j).c_str()), Form("plots/run%s/detvar/comparisons/cuts/%s/%s.pdf", _util.run_period, _util.cut_dirs.at(i).c_str(), _util.vec_hist_name.at(j).c_str()),
                            _util.cut_dirs.at(i), _util.vec_axis_label.at(j).c_str());

            }


            /*

            // Space Charge Corrected X position comparision plot
            SysVariations("h_reco_vtx_x_sce", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_vtx_x_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex X [cm]", _util.cut_dirs.at(i).c_str(), "h_reco_vtx_x_sce", file_sys_var );
    
            // Space Charge Corrected Y position comparision plot
            SysVariations("h_reco_vtx_y_sce", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_vtx_y_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Y [cm]",  _util.cut_dirs.at(i).c_str(), "h_reco_vtx_y_sce", file_sys_var);

            // Space Charge Corrected X position comparision plot
            SysVariations("h_reco_vtx_z_sce", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_vtx_z_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Z [cm]", _util.cut_dirs.at(i).c_str(), "h_reco_vtx_z_sce", file_sys_var);
        
            // Flash Time
            SysVariations("h_reco_flash_time", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_flash_time.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Flash Time [#mus]", _util.cut_dirs.at(i).c_str(), "h_reco_flash_time", file_sys_var);

            // Leading Shower Phi
            SysVariations("h_reco_leading_shower_phi", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_leading_shower_phi.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Phi [degrees]", _util.cut_dirs.at(i).c_str(), "h_reco_leading_shower_phi", file_sys_var);

            // Leading Shower Theta
            SysVariations("h_reco_leading_shower_theta", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_leading_shower_theta.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Theta [degrees]", _util.cut_dirs.at(i).c_str(), "h_reco_leading_shower_theta", file_sys_var);

            // Shower Multiplicty
            SysVariations("h_reco_shower_multiplicity", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_shower_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Multiplicty", _util.cut_dirs.at(i).c_str(), "h_reco_shower_multiplicity", file_sys_var);

            // Track Multiplicty
            SysVariations("h_reco_track_multiplicity", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_track_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Track Multiplicity", _util.cut_dirs.at(i).c_str(), "h_reco_track_multiplicity", file_sys_var);

            // Topological Score
            SysVariations("h_reco_topological_score", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_topological_score.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Topological Score", _util.cut_dirs.at(i).c_str(), "h_reco_topological_score", file_sys_var);
            
            // Shower Track Fitter dedx Y plane
            SysVariations("h_reco_shr_tkfit_dedx_y", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_shr_tkfit_dedx_y.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Collection Plane dEdx (track fitter) [MeV/cm]", _util.cut_dirs.at(i).c_str(), "h_reco_shr_tkfit_dedx_y", file_sys_var);

            // Reco Electron Neutrino E
            SysVariations("h_reco_nu_e", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_nu_e.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Neutrino Energy [GeV]", _util.cut_dirs.at(i).c_str(), "h_reco_nu_e", file_sys_var);

            // Leading Shower Energy
            SysVariations("h_reco_shower_energy_tot_cali", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_shower_energy_tot_cali.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Leading Shower Energy [GeV]", _util.cut_dirs.at(i).c_str(), "h_reco_shower_energy_tot_cali", file_sys_var);
    
            // --------------------- new variables asked by Krish

            // Software Trigger
            SysVariations("h_reco_softwaretrig", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_softwaretrigger.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Software Trigger", _util.cut_dirs.at(i).c_str(), "h_reco_softwaretrig", file_sys_var);

            // nslice
            SysVariations("h_reco_nslice", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_slice.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "nslice", _util.cut_dirs.at(i).c_str(), "h_reco_nslice", file_sys_var);

            // Shower Score
            SysVariations("h_reco_shower_score", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_shower_score.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Score", _util.cut_dirs.at(i).c_str(), "h_reco_shower_score", file_sys_var);

            // dEdx max
            SysVariations("h_reco_shr_tkfit_dedx_max", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_dedx_max.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "dE/dx max", _util.cut_dirs.at(i).c_str(), "h_reco_shr_tkfit_dedx_max", file_sys_var);
    
            // dEdx max with tracks
            SysVariations("h_reco_shr_tkfit_dedx_max_with_tracks", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_dedx_max_with_track.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "dE/dx max with tracks", _util.cut_dirs.at(i).c_str(), "h_reco_shr_tkfit_dedx_max_with_tracks", file_sys_var);

            // dEdx max no tracks
            SysVariations("h_reco_shr_tkfit_dedx_max_no_tracks", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_dedx_max_no_tracks.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "dE/dx max no tracks", _util.cut_dirs.at(i).c_str(), "h_reco_shr_tkfit_dedx_max_no_tracks", file_sys_var);

            // Shower Distance
            SysVariations("h_reco_shower_to_vtx_dist", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_shower_to_vtx_dist.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower to Vertex Distance [cm]", _util.cut_dirs.at(i).c_str(), "h_reco_shower_to_vtx_dist", file_sys_var);

            // Hits Ratio
            SysVariations("h_reco_hits_ratio", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_hits_ratio.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Hits Ratio", _util.cut_dirs.at(i).c_str(), "h_reco_hits_ratio", file_sys_var);

            // Cosmic IPAll 3D
            SysVariations("h_reco_CosmicIPAll3D", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_CosmicIPAll3D.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Cosmic IPall 3D", _util.cut_dirs.at(i).c_str(), "h_reco_CosmicIPAll3D", file_sys_var);

            // Contained Fraction
            SysVariations("h_reco_contained_fraction", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_contained_fraction.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Contained Fraction", _util.cut_dirs.at(i).c_str(), "h_reco_contained_fraction", file_sys_var);

            // Shower Molier Average
            SysVariations("h_reco_shrmoliereavg", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_shrmoliereavg.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Moliere Average", _util.cut_dirs.at(i).c_str(), "h_reco_shrmoliereavg", file_sys_var);

            // Shower Energy Cali
            SysVariations("h_reco_shower_energy_cali", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_shower_energy_cali.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Energy Cali [GeV]", _util.cut_dirs.at(i).c_str(), "h_reco_shower_energy_cali", file_sys_var);

            // Flash PE
            SysVariations("h_reco_flash_pe", Form("plots/run%s/systvar/comparisons/cuts/%s/reco_flash_pe.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Flash PE", _util.cut_dirs.at(i).c_str(), "h_reco_flash_pe", file_sys_var);
				
                        
            // Create the directory
            _util.CreateDirectory("/detvar/comparisons/cuts/" + _util.cut_dirs.at(i));

	        std::cout << "Printing detvar" << std::endl;

            // Space Charge Corrected X position comparision plot
            PlotVariations("h_reco_vtx_x_sce", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_vtx_x_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex X [cm]");

            // Space Charge Corrected Y position comparision plot
            PlotVariations("h_reco_vtx_y_sce", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_vtx_y_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Y [cm]");


            // Space Charge Corrected X position comparision plot
            PlotVariations("h_reco_vtx_z_sce", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_vtx_z_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Z [cm]");
        
            // Flash Time
            PlotVariations("h_reco_flash_time", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_flash_time.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Flash Time [#mus]");

            // Leading Shower Phi
            PlotVariations("h_reco_leading_shower_phi", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_leading_shower_phi.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Phi [degrees]");

            // Leading Shower Theta
            PlotVariations("h_reco_leading_shower_theta", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_leading_shower_theta.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Theta [degrees]");


            // Shower Multiplicty
            PlotVariations("h_reco_shower_multiplicity", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_shower_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Multiplicty");

            // Track Multiplicty
            PlotVariations("h_reco_track_multiplicity", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_track_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Track Multiplicity");

            
            // Topological Score
            PlotVariations("h_reco_topological_score", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_topological_score.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Topological Score");

            
            // Shower Track Fitter dedx Y plane
            PlotVariations("h_reco_shr_tkfit_dedx_y", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_shr_tkfit_dedx_y.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Collection Plane dEdx (track fitter) [MeV/cm]");

            // Reco Electron Neutrino E
            PlotVariations("h_reco_nu_e", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_nu_e.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Neutrino Energy [GeV]");

            // Leading Shower Energy
            PlotVariations("h_reco_shower_energy_tot_cali", Form("plots/run%s/detvar/comparisons/cuts/%s/reco_shower_energy_tot_cali.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Leading Shower Energy [GeV]");

                            */

        }
        // Ext mode
        else if (mode == "ext"){

            // Create the directory
            _util.CreateDirectory("/ext/comparisons/cuts/" + _util.cut_dirs.at(i));

            // Space Charge Corrected X position comparision plot
            PlotVariationsEXT("h_reco_vtx_x_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_x_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex X [cm]");

            // Space Charge Corrected Y position comparision plot
            PlotVariationsEXT("h_reco_vtx_y_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_y_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Y [cm]");


            // Space Charge Corrected X position comparision plot
            PlotVariationsEXT("h_reco_vtx_z_sce", Form("plots/run%s/ext/comparisons/cuts/%s/reco_vtx_z_sce.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reco Vertex Z [cm]");
        
            // Leading Shower Phi
            PlotVariationsEXT("h_reco_leading_shower_phi", Form("plots/run%s/ext/comparisons/cuts/%s/reco_leading_shower_phi.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Phi [degrees]");

            // Leading Shower Theta
            PlotVariationsEXT("h_reco_leading_shower_theta", Form("plots/run%s/ext/comparisons/cuts/%s/reco_leading_shower_theta.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Leading Shower Theta [degrees]");


            // Shower Multiplicty
            PlotVariationsEXT("h_reco_shower_multiplicity", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shower_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Shower Multiplicty");

            // Track Multiplicty
            PlotVariationsEXT("h_reco_track_multiplicity", Form("plots/run%s/ext/comparisons/cuts/%s/reco_track_multiplicity.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Track Multiplicity");

            
            // Topological Score
            PlotVariationsEXT("h_reco_topological_score", Form("plots/run%s/ext/comparisons/cuts/%s/reco_topological_score.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Topological Score");

            
            // Shower Track Fitter dedx Y plane
            PlotVariationsEXT("h_reco_shr_tkfit_dedx_y", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shr_tkfit_dedx_y.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Collection Plane dEdx (track fitter) [MeV/cm]");

            // Reco Electron Neutrino E
            PlotVariationsEXT("h_reco_nu_e", Form("plots/run%s/ext/comparisons/cuts/%s/reco_nu_e.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Neutrino Energy [GeV]");

            // Leading Shower Energy
            PlotVariationsEXT("h_reco_shower_energy_tot_cali", Form("plots/run%s/ext/comparisons/cuts/%s/reco_shower_energy_tot_cali.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Reconstructed Leading Shower Energy [GeV]");

            // Flash PE
            PlotVariationsEXT("h_reco_flash_pe", Form("plots/run%s/ext/comparisons/cuts/%s/reco_flash_pe.pdf", _util.run_period, _util.cut_dirs.at(i).c_str()),
                            _util.cut_dirs.at(i), "Flash PE [PE]");

        }
        // dont know what your trying to configure tbh ;)
        else {
            std::cout << "something is going wrong here..." << std::endl;
        }
        
    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotVariations(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name ){

    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    hist.resize(k_vars_MAX);
    hist_ratio.resize(k_vars_MAX);

    TCanvas * c      = new TCanvas();
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    _util.SetTPadOptions(topPad, bottomPad );

    // Loop over the variations and get the histograms
    for (unsigned int k=0; k < f_vars.size(); k++){
        
        // Loop over the classifications and get the histograms
        for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

            // Only want the MC piece -- may want to add in dirt too? -- will need to separately scale that hitogram though
            if ( i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt ) continue;

            // Get all the MC histograms and add them together
            TH1D *h_temp;

            _util.GetHist(f_vars.at(k), h_temp, Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            
            // First case so clone the histogram
            if (i == 0) hist.at(k) = (TH1D*) h_temp->Clone("h_sum_hist");
            else hist.at(k)->Add(h_temp, 1);
        }
        
    }

    // Now scale the histograms to POT
    for (unsigned int y=0; y < hist.size(); y++ ){
        double scale_fact = POT_v.at(k_CV) / POT_v.at(y);
        // std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        if (y == k_CV){
            // Clone a histogram to plot the CV error as a grey band
            h_error_hist = (TH1D*) hist.at(k_CV)->Clone("h_error_hist");
            h_error_hist->SetFillColorAlpha(12, 0.15);
            h_error_hist->SetLineWidth(2);
            h_error_hist->SetLineColor(kBlack);
            h_error_hist->GetYaxis()->SetTitle("Entries");
            h_error_hist->GetYaxis()->SetTitleFont(46);
            h_error_hist->GetYaxis()->SetTitleSize(13);
            h_error_hist->Draw("E2");
        }


        // Save clones of the histograms for doing the ratios
        hist_ratio.at(y) = (TH1D*) hist.at(y)->Clone(Form("h_ratio_%s", var_string.at(y).c_str()));
        hist_ratio.at(y)->Divide(hist.at(k_CV)); 

        // Set the customisation of the histogram
        SetVariationProperties(hist.at(y), y);

        // Draw the histograms
        if (y == k_CV) hist.at(y)->Draw("hist, same");
        else hist.at(y)->Draw("hist,E, same");
    }

    // Legend
    TLegend *leg = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_error_hist, "CV",   "lf");
    leg->AddEntry(hist.at(k_bnb_diffusion), "BNB Diffusion", "l");
    leg->Draw();

    // Now draw the ratio
    bottomPad->cd();

    for (unsigned int k =0; k < hist_ratio.size(); k++){

        SetVariationProperties(hist_ratio.at(k), k);
        hist_ratio.at(k)->SetLineWidth(1);


        hist_ratio.at(k)->GetXaxis()->SetLabelSize(12);
        hist_ratio.at(k)->GetXaxis()->SetLabelFont(43); 
        hist_ratio.at(k)->GetYaxis()->SetLabelSize(11);
        hist_ratio.at(k)->GetYaxis()->SetLabelFont(43);
        hist_ratio.at(k)->GetXaxis()->SetTitleOffset(3.0);
        hist_ratio.at(k)->GetXaxis()->SetTitleSize(17);
        hist_ratio.at(k)->GetXaxis()->SetTitleFont(46);
        hist_ratio.at(k)->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        hist_ratio.at(k)->GetYaxis()->SetRangeUser(0., 2.0);
        hist_ratio.at(k)->GetYaxis()->SetTitle("Variation / CV");
        hist_ratio.at(k)->GetYaxis()->SetTitleSize(13);
        hist_ratio.at(k)->GetYaxis()->SetTitleFont(44);
        hist_ratio.at(k)->GetYaxis()->SetTitleOffset(1.5);
        hist_ratio.at(k)->SetTitle(" ");
        hist_ratio.at(k)->GetXaxis()->SetTitle(x_axis_name);
        
        if (k == 0) hist_ratio.at(k)->Draw("hist,same");
        else hist_ratio.at(k)->Draw("hist,E,same");
    }

    // Draw the error hist 
    TH1D* h_ratio_error = (TH1D*) h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");


    // Draw the run period on the plot
    // _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);
    
    c->Print(print_name);



}
// -----------------------------------------------------------------------------
void SystematicsHelper::SysVariations(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name, std::string folder_name, std::string plot_name, TFile *root_output){

    // last updated on Sept 18 by Marina Reggiani-Guzzo

    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_diff; // The vector of histogram differentes between CV and the vatiation (variation-CV)
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    hist.resize(k_vars_MAX);
    hist_diff.resize(k_vars_MAX);
    hist_ratio.resize(k_vars_MAX);

    TCanvas * c      = new TCanvas();
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    _util.SetTPadOptions(topPad, bottomPad );

    // Loop over the variations and get the histograms
    for (unsigned int k=0; k < f_vars.size(); k++){
        
        // Loop over the classifications and get the histograms
        for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

            // Only want the MC piece -- may want to add in dirt too? -- will need to separately scale that hitogram though
            if ( i == _util.k_leg_data || i == _util.k_leg_ext || i == _util.k_leg_dirt ) continue;

            // Get all the MC histograms and add them together
            TH1D *h_temp;

            _util.GetHist(f_vars.at(k), h_temp, Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            
            // First case so clone the histogram
            if (i == 0) hist.at(k) = (TH1D*) h_temp->Clone("h_sum_hist");
            else hist.at(k)->Add(h_temp, 1);
        }
        
    }

    // Now scale the histograms to POT
    for (unsigned int y=0; y < hist.size(); y++ ){
        double scale_fact = POT_v.at(k_CV) / POT_v.at(y);
        // std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        if (y == k_CV){
            // Clone a histogram to plot the CV error as a grey band
            h_error_hist = (TH1D*) hist.at(k_CV)->Clone("h_error_hist");
            h_error_hist->SetFillColorAlpha(12, 0.15);
            h_error_hist->SetLineWidth(2);
            h_error_hist->SetLineColor(kBlack);
            h_error_hist->GetYaxis()->SetTitle("Entries");
            h_error_hist->GetYaxis()->SetTitleFont(46);
            h_error_hist->GetYaxis()->SetTitleSize(13);
            h_error_hist->Draw("E2");
        }

        // calculate difference between cv and the histogram with variation
        hist_diff.at(y) = (TH1D*) hist.at(y)->Clone(Form("h_diff_%s", var_string.at(y).c_str()));
        hist_diff.at(y)->Add(hist.at(k_CV),-1); 

        // Save clones of the histograms for doing the ratios
        hist_ratio.at(y) = (TH1D*) hist_diff.at(y)->Clone(Form("h_ratio_%s", var_string.at(y).c_str()));
        hist_ratio.at(y)->Divide(hist.at(k_CV));

        // Set the customisation of the histogram
        SetVariationProperties(hist.at(y), y);

        // Draw the histograms
        if (y == k_CV) hist.at(y)->Draw("hist, same");
        else hist.at(y)->Draw("hist,E, same");
    }

    // Legend
    TLegend *leg = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_error_hist, "CV",   "lf");
    leg->AddEntry(hist.at(k_bnb_diffusion), "BNB Diffusion", "l");
    leg->Draw();

    // Now draw the ratio
    bottomPad->cd();

    for (unsigned int k =0; k < hist_ratio.size(); k++){

        SetVariationProperties(hist_ratio.at(k), k);
        hist_ratio.at(k)->SetLineWidth(1);


        hist_ratio.at(k)->GetXaxis()->SetLabelSize(12);
        hist_ratio.at(k)->GetXaxis()->SetLabelFont(43); 
        hist_ratio.at(k)->GetYaxis()->SetLabelSize(11);
        hist_ratio.at(k)->GetYaxis()->SetLabelFont(43);
        hist_ratio.at(k)->GetXaxis()->SetTitleOffset(3.0);
        hist_ratio.at(k)->GetXaxis()->SetTitleSize(17);
        hist_ratio.at(k)->GetXaxis()->SetTitleFont(46);
        hist_ratio.at(k)->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        hist_ratio.at(k)->GetYaxis()->SetRangeUser(-3.0, 3.0);
        hist_ratio.at(k)->GetYaxis()->SetTitle("(Variation-CV) / CV");
        hist_ratio.at(k)->GetYaxis()->SetTitleSize(13);
        hist_ratio.at(k)->GetYaxis()->SetTitleFont(44);
        hist_ratio.at(k)->GetYaxis()->SetTitleOffset(1.5);
        hist_ratio.at(k)->SetTitle(" ");
        hist_ratio.at(k)->GetXaxis()->SetTitle(x_axis_name);
        
        if (k == 0) hist_ratio.at(k)->Draw("hist,same");
        else {

		hist_ratio.at(k)->Draw("hist,E,same");

	}

   }

	// loop over the variations and save them into different directories
	// this root file is going to be used later on to calculate the sys uncertainty	
	for (unsigned int k=0; k < f_vars.size(); k++){
		
		if(!root_output->GetDirectory(Form("%s/%s", folder_name.c_str(), var_string.at(k).c_str()))) {
			root_output->mkdir(Form("%s/%s", folder_name.c_str(), var_string.at(k).c_str())); // if the directory does not exist, create it
		}

		root_output->cd(Form("%s/%s", folder_name.c_str(), var_string.at(k).c_str())); // open the directory
	
		hist_ratio.at(k)->SetDirectory(gDirectory); // set in which dir the hist_ratio.at(k) is going to be written
		hist_ratio.at(k)->Write(Form("%s", plot_name.c_str()), TObject::kOverwrite);  // write the histogram to the file
	
	}

    

    // Draw the error hist 
    /*
    TH1D* h_ratio_error = (TH1D*) h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");
    */


    // Draw the run period on the plot
    // _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Add the weight labels
    // Draw_WeightLabels(c);
    
    c->Print(print_name);

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotVariationsEXT(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name ){

    gStyle->SetOptStat(0);

    std::vector<TH1D*> hist; // The vector of histograms from the file for the plot
    std::vector<TH1D*> hist_ratio; // The vector of histogram ratios from the file for the plot
    TH1D * h_error_hist;
    
    hist.resize(var_string.size());
    hist_ratio.resize(var_string.size());

    TCanvas * c      = new TCanvas();
    TPad * topPad    = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    TPad * bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);

    _util.SetTPadOptions(topPad, bottomPad );


    // Loop over the variations and get the histograms
    for (unsigned int k=0; k < f_vars.size(); k++){
        
        // Loop over the classifications and get the histograms
        for (unsigned int i=0; i <_util.classification_dirs.size(); i++){

            // Only want the EXT piece
            if ( i != _util.k_leg_ext ) continue;

            _util.GetHist(f_vars.at(k), hist.at(k), Form("Stack/%s/%s/%s_%s_%s", cut_name.c_str(), _util.classification_dirs.at(i).c_str(), hist_name.c_str(), cut_name.c_str(), _util.classification_dirs.at(i).c_str()));
            
        }
        
    }


    // Now scale the histograms to POT
    for (unsigned int y=0; y < hist.size(); y++ ){
        double scale_fact = POT_v.at(k_NuMI) / POT_v.at(y);
        // std::cout << "scale factor: " << scale_fact << std::endl;
        hist.at(y)->Scale(scale_fact);

        // Do area norm
        std::cout << "BNB to NuMI scale Factor: " << hist.at(k_NuMI)->Integral() / hist.at(y)->Integral() << std::endl;
        hist.at(y)->Scale(hist.at(k_NuMI)->Integral() / hist.at(y)->Integral() );
        
        if (y == k_NuMI){
            // Clone a histogram to plot the CV error as a grey band
            h_error_hist = (TH1D*) hist.at(k_NuMI)->Clone("h_error_hist");
            h_error_hist->SetFillColorAlpha(12, 0.15);
            h_error_hist->SetLineWidth(2);
            h_error_hist->SetLineColor(kBlack);
            h_error_hist->GetYaxis()->SetTitle("Entries");
            h_error_hist->GetYaxis()->SetTitleFont(46);
            h_error_hist->GetYaxis()->SetTitleSize(13);
            h_error_hist->Draw("E2");
        }


        // Save clones of the histograms for doing the ratios
        hist_ratio.at(y) = (TH1D*) hist.at(y)->Clone(Form("h_ratio_%s", var_string.at(y).c_str()));
        hist_ratio.at(y)->Divide(hist.at(k_NuMI)); 

        // Set the customisation of the histogram
        SetVariationProperties(hist.at(y), y);

        // Draw the histograms
        if (y == k_NuMI) hist.at(y)->Draw("hist, same");
        else hist.at(y)->Draw("hist,E, same");

    }

    // Legend
    TLegend *leg = new TLegend(0.8, 0.91, 0.95, 0.32);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_error_hist, "NuMI EXT",   "lf");
    leg->AddEntry(hist.at(k_BNB), "BNB EXT", "l");
    leg->Draw();

    // Now draw the ratio
    bottomPad->cd();

    for (unsigned int k =0; k < hist_ratio.size(); k++){

        SetVariationProperties(hist_ratio.at(k), k);
        hist_ratio.at(k)->SetLineWidth(1);


        hist_ratio.at(k)->GetXaxis()->SetLabelSize(12);
        hist_ratio.at(k)->GetXaxis()->SetLabelFont(43); 
        hist_ratio.at(k)->GetYaxis()->SetLabelSize(11);
        hist_ratio.at(k)->GetYaxis()->SetLabelFont(43);
        hist_ratio.at(k)->GetXaxis()->SetTitleOffset(3.0);
        hist_ratio.at(k)->GetXaxis()->SetTitleSize(17);
        hist_ratio.at(k)->GetXaxis()->SetTitleFont(46);
        hist_ratio.at(k)->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        hist_ratio.at(k)->GetYaxis()->SetRangeUser(0., 2.0);
        hist_ratio.at(k)->GetYaxis()->SetTitle("NuMI / BNB");
        hist_ratio.at(k)->GetYaxis()->SetTitleSize(13);
        hist_ratio.at(k)->GetYaxis()->SetTitleFont(44);
        hist_ratio.at(k)->GetYaxis()->SetTitleOffset(1.5);
        hist_ratio.at(k)->SetTitle(" ");
        hist_ratio.at(k)->GetXaxis()->SetTitle(x_axis_name);
        
        if (k == 0) hist_ratio.at(k)->Draw("hist,same");
        else hist_ratio.at(k)->Draw("hist,E,same");
    }

    // Draw the error hist 
    TH1D* h_ratio_error = (TH1D*) h_error_hist->Clone("h_ratio_error");
    h_ratio_error->Divide(h_ratio_error);
    h_ratio_error->Draw("e2, same");


    // Draw the run period on the plot
    _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

    // Draw area normalisation
    Draw_Area_Norm(c);
    
    c->Print(print_name);



}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseReweightingMode(){

    gStyle->SetOptStat(0);

    // Load in the input file
    // Should we add more protection to this command??
    f_nuexsec = new TFile( Form("files/crosssec_run%s.root", _util.run_period ), "READ");

    InitialsePlotCV();

    // Now lets initialse the vectors which will store the total uncertainties
    InitialiseUncertaintyVectors();


    // Loop over the cross-section variables
    for (unsigned int var = 0; var <  vars.size(); var++){

        // Comparison plots for data to MC
        CompareVariationXSec("RPA",              var, "RPA" );
        CompareVariationXSec("CCMEC",            var, "CC MEC" );
        CompareVariationXSec("AxFFCCQE",         var, "Ax FF CCQE" );
        CompareVariationXSec("VecFFCCQE",        var, "Vec FF CCQE" );
        CompareVariationXSec("DecayAngMEC",      var, "Decay Ang MEC" );
        CompareVariationXSec("ThetaDelta2Npi",   var, "Theta Delta 2N #pi" );
        CompareVariationXSec("ThetaDelta2NRad",  var, "Theta Delta 2N Rad" );
        CompareVariationXSec("RPA_CCQE_Reduced", var, "RPA CCQE Reduced" );
        CompareVariationXSec("NormCCCOH",        var, "Norm CC COH" );
        CompareVariationXSec("NormNCCOH",        var, "Norm NC COH" );

        // Plot the interaction unisims
        PlotReweightingModeUnisim("RPA",              var, "RPA" );
        PlotReweightingModeUnisim("CCMEC",            var, "CC MEC" );
        PlotReweightingModeUnisim("AxFFCCQE",         var, "Ax FF CCQE" );
        PlotReweightingModeUnisim("VecFFCCQE",        var, "Vec FF CCQE" );
        PlotReweightingModeUnisim("DecayAngMEC",      var, "Decay Ang MEC" );
        PlotReweightingModeUnisim("ThetaDelta2Npi",   var, "Theta Delta 2N #pi" );
        PlotReweightingModeUnisim("ThetaDelta2NRad",  var, "Theta Delta 2N Rad" );
        PlotReweightingModeUnisim("RPA_CCQE_Reduced", var, "RPA CCQE Reduced" );
        PlotReweightingModeUnisim("NormCCCOH",        var, "Norm CC COH" );
        PlotReweightingModeUnisim("NormNCCOH",        var, "Norm NC COH" );

        // Dirt
        PlotReweightingModeUnisim("Dirt",        var, "Dirt" );

        // Plot the beamline unisims
        PlotReweightingModeUnisim("Horn_curr",          var, "Horn Current" );
        PlotReweightingModeUnisim("Horn1_x",            var, "Horn 1 x" );
        PlotReweightingModeUnisim("Horn1_y",            var, "Horn 1 y" );
        PlotReweightingModeUnisim("Beam_spot",          var, "Beam Spot Size" );
        PlotReweightingModeUnisim("Horn2_x",            var, "Horn 2 x" );
        PlotReweightingModeUnisim("Horn2_y",            var, "Horn 2 y" );
        PlotReweightingModeUnisim("Horn_Water",         var, "Horns Water" );
        PlotReweightingModeUnisim("Beam_shift_x",       var, "Beam shift x" );
        PlotReweightingModeUnisim("Beam_shift_y",       var, "Beam shift y" );
        PlotReweightingModeUnisim("Target_z",           var, "Target z" );
        PlotReweightingModeUnisim("Horn1_refined_descr",var, "Horn 1 Refined Desc." );
        PlotReweightingModeUnisim("Decay_pipe_Bfield",  var, "Decay pipe Bfield" );
        PlotReweightingModeUnisim("Old_Horn_Geometry",  var, "Old Horn Geometry" );

        // Plot the multisims
        PlotReweightingModeMultisim("weightsGenie", var,  "GENIE", 600);
        PlotReweightingModeMultisim("weightsReint", var,  "Geant Reinteractions", 1000);
        PlotReweightingModeMultisim("weightsPPFX",  var,  "PPFX", 600);
        
    }

    // Get the statistical uncertainties
    FillStatVector();

    // Add 2% POT counting uncertainty to the tot sys error vector
    FillPOTCountingVector();

    // Compare the MC and Data X-Section
    CompareCVXSec();
    CompareCVXSecNoRatio();

    PrintUncertaintySummary();

}
// -----------------------------------------------------------------------------
void SystematicsHelper::SetLabelName(std::string label, std::string &label_up, std::string &label_dn){

    if       (label == "Horn_curr"          ){
        label_up = "Horn_p2kA";
        label_dn = "Horn_m2kA";
    }
    else if  (label == "Horn1_x"       ){
        label_up = label + "_p3mm";
        label_dn = "Horm1_x_m3mm";
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
        label_dn = "Horm2_x_m3mm";
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
    else if  (label == "Horn1_refined_descr" || label == "Decay_pipe_Bfield" || label == "Old_Horn_Geometry"){
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
void SystematicsHelper::PlotReweightingModeUnisim(std::string label, int var, std::string label_pretty){

    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe;
    
    // Resize to the number of universes
    h_universe.resize(2);
    h_universe.at(k_up).resize(xsec_types.size());
    h_universe.at(k_dn).resize(xsec_types.size());

    // Now get the histograms
    std::string label_up = label + "up";
    std::string label_dn = label + "dn";

    // Set the Unisim up down variation name
    SetLabelName(label, label_up, label_dn);

    // Check if its just a single on/off type variation
    // This case we dont want to plot the up/dn, but just once
    bool single_var = false;
    if (label_up == label_dn) single_var = true;

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        _util.GetHist(f_nuexsec, h_universe.at(k_up).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), _util.run_period, label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, h_universe.at(k_dn).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), _util.run_period, label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
    }

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c;
    
    // Now we want to draw them
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        c = new TCanvas("c", "c", 500, 500);
        topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        _util.SetTPadOptions(topPad, bottomPad);
        // topPad->SetRightMargin(0.10 );
        // bottomPad->SetRightMargin(0.10 );
        
        h_universe.at(k_up).at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));
        h_universe.at(k_up).at(k)->GetXaxis()->SetTitle("");
        h_universe.at(k_up).at(k)->GetXaxis()->SetLabelSize(0);

        h_universe.at(k_up).at(k)->Draw("hist");
        cv_hist_vec.at(var).at(k)->Draw("hist,same");
        if (!single_var) h_universe.at(k_dn).at(k)->Draw("hist,same");

        c->Update();

        double scale_val = h_universe.at(k_up).at(k)->GetMaximum();
        if (scale_val < cv_hist_vec.at(var).at(k)->GetMaximum()) scale_val = cv_hist_vec.at(var).at(k)->GetMaximum();
        if (scale_val <  h_universe.at(k_dn).at(k)->GetMaximum()) scale_val =  h_universe.at(k_dn).at(k)->GetMaximum();

        h_universe.at(k_up).at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        // FIxed scaling for differential cross section
        if (vars.at(var) != "integrated"){
            // h_universe.at(k_up).at(k)->GetYaxis()->SetRangeUser(0, 0.5e-39);
        }

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

        bottomPad->cd();
        
        // Up ratio to CV
        TH1D* h_err_up = (TH1D *)h_universe.at(k_up).at(k)->Clone("h_ratio_up");
        h_err_up->Add(cv_hist_vec.at(var).at(k), -1);
        h_err_up->Divide(cv_hist_vec.at(var).at(k));
        h_err_up->SetLineWidth(2);
        h_err_up->SetLineColor(kGreen+2);
        h_err_up->Scale(100);

        h_err_up->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
        
        // Down ratio to CV
        TH1D* h_err_dn = (TH1D *)h_universe.at(k_dn).at(k)->Clone("h_ratio_dn");
        h_err_dn->Add(cv_hist_vec.at(var).at(k), -1);
        h_err_dn->Divide(cv_hist_vec.at(var).at(k));
        h_err_dn->SetLineWidth(2);
        h_err_dn->SetLineColor(kRed+2);
        h_err_dn->Scale(100);

        TH1D* h_err = (TH1D *)cv_hist_vec.at(var).at(k)->Clone("h_ratio");
        h_err->Divide(cv_hist_vec.at(var).at(k));

        SetRatioOptions(h_err_up);
        h_err_up->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
        h_err_up->GetYaxis()->SetTitle("\% change from CV");
        h_err_up->Draw("hist,same");

        if (!single_var) h_err_dn->Draw("hist,same");

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
        
        h_err->Draw("hist,same");

        FillSysVector(label, var, k, h_err_up, h_err_dn);

        c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(), _util.run_period, label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;
        delete h_err;
        delete h_err_up;
        delete h_err_dn;
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::PlotReweightingModeMultisim(std::string label, int var, std::string label_pretty, int universes){

    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe; // Universe, <gen/sig/xsec etc>
    std::vector<std::vector<TH1D*>> h_err;


    // Clone the CV histograms so we can work with them without changing them
    std::vector<TH1D*> cv_hist_vec_clone;
    cv_hist_vec_clone.resize(xsec_types.size());
    
    for (unsigned int h_index = 0; h_index < cv_hist_vec_clone.size(); h_index++){
        cv_hist_vec_clone.at(h_index) = (TH1D*)cv_hist_vec.at(var).at(h_index)->Clone(Form("h_%s_clone", xsec_types.at(h_index).c_str() ));
        // Customise
        cv_hist_vec_clone.at(h_index)->SetLineWidth(2);
        cv_hist_vec_clone.at(h_index)->SetLineColor(kBlack);
    }

    
    // Resize to the number of universes
    h_universe.resize(universes);
    
    // Resize each universe
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        h_universe.at(uni).resize(xsec_types.size());
    }
    
    
    // Get the histograms and customise a bit
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            _util.GetHist(f_nuexsec, h_universe.at(uni).at(k), Form( "%s/%s/h_run%s_%s_%i_%s_%s", label.c_str(), vars.at(var).c_str(), _util.run_period, label.c_str(), uni ,vars.at(var).c_str(), xsec_types.at(k).c_str()));

            // Customise
            h_universe.at(uni).at(k)->SetLineWidth(1);
            h_universe.at(uni).at(k)->SetLineColor(kAzure+5);
            
        }
    }

    // Get the covariance matrix
    CalcCovariance(label, var, h_universe );

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



    TCanvas *c;
    TPad *topPad;
    TPad *bottomPad;

    // Now we want to draw them
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){

        c = new TCanvas("c", "c", 500, 500);
        topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
        bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
        _util.SetTPadOptions(topPad, bottomPad);

        h_universe.at(0).at(k)->SetTitle(Form("%s", xsec_types_pretty.at(k).c_str() ));

        double scale_val = h_universe.at(0).at(k)->GetMaximum();
        
        if (scale_val <  h_universe.at(k_dn).at(k)->GetMaximum()) scale_val =  h_universe.at(k_dn).at(k)->GetMaximum();

        // Loop over universes
        for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
            h_universe.at(uni).at(k)->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            h_universe.at(uni).at(k)->Draw("hist,same");
            if (scale_val < h_universe.at(uni).at(k)->GetMaximum()) scale_val = h_universe.at(uni).at(k)->GetMaximum();
            h_universe.at(uni).at(k)->GetXaxis()->SetTitle("");
            h_universe.at(uni).at(k)->GetXaxis()->SetLabelSize(0);

        }

        if (xsec_types.at(k) != "data_xsec") {
            cv_hist_vec_clone.at(k)->SetLineColor(kBlack);
            cv_hist_vec_clone.at(k)->SetFillStyle(0);
            cv_hist_vec_clone.at(k)->SetMarkerColor(kBlack);
            cv_hist_vec_clone.at(k)->Draw("E2,hist,same");
        }
        else { 
            cv_hist_vec_clone.at(k)->Draw("E,same");
        }

        h_universe.at(0).at(k)->GetYaxis()->SetRangeUser(0, scale_val*1.2);

        if (k == k_xsec_mcxsec || k == k_xsec_dataxsec){
            // h_universe.at(0).at(k)->GetYaxis()->SetRangeUser(0, 0.5e-39);
        }

        TLegend *leg = new TLegend(0.5, 0.65, 0.85, 0.8);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_universe.at(0).at(k), Form("%s", label_pretty.c_str()), "l");
        if (xsec_types.at(k) != "data_xsec") leg->AddEntry(cv_hist_vec_clone.at(k),           "CV (Sys Only)", "f");
        else                                 leg->AddEntry(cv_hist_vec_clone.at(k),           "CV (Sys Only)", "le");
        leg->Draw();

        bottomPad->cd();

        // Up percent diff to CV
        TH1D* h_err = (TH1D *)cv_hist_vec_clone.at(k)->Clone("h_err");
        
        // Loop over the bins in the up error, and set the bin content to be the percent difference
        for (int g = 1; g < h_err->GetNbinsX()+1; g++){
            h_err->SetBinContent(g, 100 * h_err->GetBinError(g)/h_err->GetBinContent(g));
        }
        h_err->SetLineWidth(2);
        h_err->SetLineColor(kAzure+5);
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
        h_err->SetTitle(" ");
        h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
        
        h_err->GetYaxis()->SetTitle("\% Uncertainty");
        h_err->SetMarkerSize(4);
        h_err->Draw("hist, text00");
        gStyle->SetPaintTextFormat("4.1f");

        FillSysVector(label, var, k, h_err, h_err);

        c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_%s.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        delete c;
        delete leg;

    }

        
    
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareCVXSec(){

    std::vector<std::string> error_type = {"stat", "sys", "tot"};

    // Loop over the variables
    for (unsigned int var = 0; var < vars.size(); var++){

        // Loop over the error labels
        for (unsigned int err_lab = 0; err_lab < error_type.size(); err_lab++){

            TH1D* h_dataxsec     = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_temp");
            TH1D* h_dataxsec_tot = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_tot_temp"); // For total uncertainty
            TH1D* h_mcxsec       = (TH1D*) cv_hist_vec.at(var).at(k_xsec_mcxsec)  ->Clone("h_mx_xsec_temp");

            TPad *topPad;
            TPad *bottomPad;

            TCanvas * c = new TCanvas("c", "c", 500, 500);
            topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
            bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
            _util.SetTPadOptions(topPad, bottomPad);

            h_dataxsec->SetLineColor(kBlack);
            h_mcxsec  ->SetLineColor(kRed+2);

            // h_dataxsec->GetYaxis()->SetRangeUser(0, 0.5e-39);
            if (vars.at(var) == "integrated") h_dataxsec->GetYaxis()->SetRangeUser(0.5e-39, 2.5e-39);
            else h_dataxsec->GetYaxis()->SetRangeUser(0.0e-39, 0.5e-39);

            h_dataxsec->GetYaxis()->SetLabelSize(0.04);
            h_dataxsec->GetYaxis()->SetTitleSize(14);
            h_dataxsec->GetYaxis()->SetTitleFont(44);
            h_dataxsec->GetYaxis()->SetTitleOffset(1.5);
            h_dataxsec->GetXaxis()->SetTitle("");
            h_dataxsec->GetXaxis()->SetLabelSize(0);
            h_dataxsec->SetMarkerStyle(20);
            h_dataxsec->SetMarkerSize(0.5);
            h_dataxsec->SetMinimum(0);
            
            // Rewrite the errors for data to sys
            if (error_type.at(err_lab) == "sys"){
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            // Overwrite error to stat + sys
            else {
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec_tot->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin) + v_stat_total.at(var).at(k_xsec_dataxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            
            
            
            h_dataxsec->Draw("E,X0");
            if (error_type.at(err_lab) == "tot"){
                h_dataxsec_tot->Draw("E1,same,X0");
                h_dataxsec->Draw("E1,same,X0");

            } 

            h_mcxsec->Draw("hist,same");
            
            TH1D* h_mcxsec_clone = (TH1D *)h_mcxsec->Clone("h_mc_clone");
            h_mcxsec_clone->SetFillColorAlpha(12, 0.15);
            h_mcxsec_clone->Draw("E2,same");
            

            TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if (error_type.at(err_lab) == "stat")      leg->AddEntry(h_dataxsec, "Data (Stat)", "le");
            else if (error_type.at(err_lab) == "sys") leg->AddEntry(h_dataxsec, "Data (Sys)", "le");
            else                                      leg->AddEntry(h_dataxsec, "Data (Stat+Sys)", "le");
            leg->AddEntry(h_mcxsec_clone,   "MC (Stat)", "lf");
            leg->Draw();


            bottomPad->cd();
                
            // The percent difference of mc wrt data
            TH1D* h_err = (TH1D *)h_dataxsec->Clone("h_ratio");
            h_err->Add(h_mcxsec, -1);
            h_err->Divide(h_dataxsec);
            h_err->Scale(100);
            h_err->GetYaxis()->SetTitle("Data - MC / Data [\%]");
            h_err->SetLineWidth(2);
            h_err->SetLineColor(kGreen+2);


            bottomPad->SetGridy(kFALSE);
            
            SetRatioOptions(h_err);
            h_err->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
            h_err->SetLineColor(kBlack);
            h_err->GetYaxis()->SetTitleSize(11);
            h_err->GetYaxis()->SetRangeUser(-100, 100);
            
            h_err->GetYaxis()->SetTitleOffset(2.5);
            h_err->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            if (vars.at(var) == "integrated")  h_err->GetXaxis()->SetLabelSize(0);
            h_err->SetMarkerSize(3.0);
            h_err->Draw("hist,text00");

            // Draw the run period on the plot
            _util.Draw_Run_Period(c, 0.86, 0.915, 0.86, 0.915);

            _util.Draw_Data_POT(c, Data_POT, 0.47, 0.915, 0.47, 0.915);

            
            c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison_%s.pdf", _util.run_period, vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), error_type.at(err_lab).c_str() ));
            delete c;
            delete h_dataxsec;
            delete h_mcxsec;
            delete h_err;

        
        }

    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareCVXSecNoRatio(){

    std::vector<std::string> error_type = {"stat", "sys", "tot"};

    // Loop over the variables
    for (unsigned int var = 0; var < vars.size(); var++){

        // Loop over the error labels
        for (unsigned int err_lab = 0; err_lab < error_type.size(); err_lab++){

            TH1D* h_dataxsec     = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_temp");
            TH1D* h_dataxsec_tot = (TH1D*) cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_data_xsec_tot_temp"); // For total uncertainty
            TH1D* h_mcxsec       = (TH1D*) cv_hist_vec.at(var).at(k_xsec_mcxsec)  ->Clone("h_mx_xsec_temp");

            TPad *topPad;
            TPad *bottomPad;

            TCanvas * c = new TCanvas("c", "c", 500, 500);

            h_dataxsec->SetLineColor(kBlack);
            h_mcxsec  ->SetLineColor(kRed+2);

            // h_dataxsec->GetYaxis()->SetRangeUser(0, 0.5e-39);
            if (vars.at(var) == "integrated") h_dataxsec->GetYaxis()->SetRangeUser(0.5e-39, 2.5e-39);
            else h_dataxsec->GetYaxis()->SetRangeUser(0.0e-39, 0.5e-39);

            _util.IncreaseLabelSize(h_dataxsec, c);
            if (vars.at(var) == "integrated")h_dataxsec->GetXaxis()->SetLabelSize(0);
            h_dataxsec->GetYaxis()->SetTitleSize(0.04);
            h_dataxsec->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
            h_dataxsec->SetMarkerStyle(20);
            h_dataxsec->SetMarkerSize(0.5);
            h_dataxsec->SetMinimum(0);
            
            // Rewrite the errors for data to sys
            if (error_type.at(err_lab) == "sys"){
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            // Overwrite error to stat + sys
            else {
                for (int bin = 0; bin < h_dataxsec->GetNbinsX(); bin++){
                    h_dataxsec_tot->SetBinError(bin+1, 0.01*std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec).at(bin) + v_stat_total.at(var).at(k_xsec_dataxsec).at(bin)) * h_dataxsec->GetBinContent(bin+1));
                }

            }
            
            
            
            h_dataxsec->Draw("E,X0");
            if (error_type.at(err_lab) == "tot"){
                h_dataxsec_tot->Draw("E1,same,X0");
                h_dataxsec->Draw("E1,same,X0");

            } 

            h_mcxsec->Draw("hist,same");
            
            TH1D* h_mcxsec_clone = (TH1D *)h_mcxsec->Clone("h_mc_clone");
            h_mcxsec_clone->SetFillColorAlpha(12, 0.15);
            h_mcxsec_clone->Draw("E2,same");
            

            TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if (error_type.at(err_lab) == "stat")     leg->AddEntry(h_dataxsec, "Data (Stat)", "le");
            else if (error_type.at(err_lab) == "sys") leg->AddEntry(h_dataxsec, "Data (Sys)", "le");
            else                                      leg->AddEntry(h_dataxsec, "Data (Stat+Sys)", "le");
            leg->AddEntry(h_mcxsec_clone,   "MC (Stat)", "lf");
            leg->Draw();

            // Draw the run period on the plot
            _util.Draw_Run_Period(c, 0.86, 0.92, 0.86, 0.92);

            _util.Draw_Data_POT(c, Data_POT, 0.52, 0.92, 0.52, 0.92);

            
            c->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_data_mc_comparison_%s_no_ratio.pdf", _util.run_period, vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), error_type.at(err_lab).c_str() ));
            delete c;
            delete h_dataxsec;
            delete h_mcxsec;

        
        }

    }


}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialsePlotCV(){

    // Get the CV histograms. These should stay constant througout the code

    cv_hist_vec.resize(vars.size());
    
    for (unsigned int var = 0; var < vars.size(); var++){
        cv_hist_vec.at(var).resize(xsec_types.size());
    }


    // Loop over the vars
    for (unsigned int var = 0; var < vars.size(); var++){
        
        // Loop over the typrs
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            _util.GetHist(f_nuexsec, cv_hist_vec.at(var).at(k), Form( "CV/%s/h_run%s_CV_0_%s_%s", vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), xsec_types.at(k).c_str()));

            if (cv_hist_vec.at(var).at(k) == NULL) std::cout << "Failed to get the histogram!" << std::endl;

            // Customise
            cv_hist_vec.at(var).at(k)->SetLineWidth(2);
            cv_hist_vec.at(var).at(k)->SetLineColor(kBlack);
        }

        // Create the CV directory and draw the CV
        _util.CreateDirectory("/Systematics/CV/" + vars.at(var) + "/");
    }
    
    TCanvas *c_cv;
    
    // Loop over the vars
    for (unsigned int var = 0; var < vars.size(); var++){
        
        // Loop over the types
        for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
            c_cv = new TCanvas("c", "c", 500, 500);
        
            cv_hist_vec.at(var).at(k)->Draw("hist");
            _util.IncreaseLabelSize(cv_hist_vec.at(var).at(k), c_cv);
        
            TLegend *leg = new TLegend(0.6, 0.8, 0.95, 0.9);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(cv_hist_vec.at(var).at(k),         "CV", "l");
            // leg->Draw();

            // Draw the run period on the plot
            _util.Draw_Run_Period(c_cv, 0.86, 0.92, 0.86, 0.92);

            c_cv->Print(Form("plots/run%s/Systematics/CV/%s/run%s_CV_%s_%s.pdf", _util.run_period, vars.at(var).c_str(), _util.run_period, vars.at(var).c_str(), xsec_types.at(k).c_str()));

            delete c_cv;
            delete leg;
        }
    }
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CompareVariationXSec(std::string label, int var, std::string label_pretty){

    
    // Create the directory
    _util.CreateDirectory("/Systematics/" + label + "/" + vars.at(var));

    std::vector<std::vector<TH1D*>> h_universe;
    
    // Resize to the number of universes
    h_universe.resize(2);
    h_universe.at(k_up).resize(xsec_types.size());
    h_universe.at(k_dn).resize(xsec_types.size());

    // Now get the histograms
    std::string label_up = label + "up";
    std::string label_dn = label + "dn";
    

    // Get the histograms and customise a bit
    for (unsigned int k = 0; k < cv_hist_vec.at(var).size(); k++){
        _util.GetHist(f_nuexsec, h_universe.at(k_up).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_up.c_str(), vars.at(var).c_str(), _util.run_period, label_up.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_up).at(k)->SetLineWidth(2);
        h_universe.at(k_up).at(k)->SetLineColor(kGreen+2);
        if (k == k_xsec_mcxsec) h_universe.at(k_up).at(k)->SetLineStyle(7);
        h_universe.at(k_up).at(k)->GetYaxis()->SetLabelSize(0.04);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleSize(14);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleFont(44);
        h_universe.at(k_up).at(k)->GetYaxis()->SetTitleOffset(1.5);
        
        _util.GetHist(f_nuexsec, h_universe.at(k_dn).at(k), Form( "%s/%s/h_run%s_%s_0_%s_%s", label_dn.c_str(), vars.at(var).c_str(), _util.run_period, label_dn.c_str(), vars.at(var).c_str(), xsec_types.at(k).c_str()));

        // Customise
        h_universe.at(k_dn).at(k)->SetLineWidth(2);
        h_universe.at(k_dn).at(k)->SetLineColor(kRed+2);
        if (k == k_xsec_mcxsec) h_universe.at(k_dn).at(k)->SetLineStyle(7);
    }

    TPad *topPad;
    TPad *bottomPad;
    TCanvas *c;
    
    // Now we want to draw them
    c = new TCanvas("c", "c", 500, 500);
    topPad = new TPad("topPad", "", 0, 0.3, 1, 1.0);
    bottomPad = new TPad("bottomPad", "", 0, 0.05, 1, 0.3);
    _util.SetTPadOptions(topPad, bottomPad);
    // topPad->SetRightMargin(0.10 );
    // bottomPad->SetRightMargin(0.10 );
    
    h_universe.at(k_up).at(k_xsec_dataxsec)->SetTitle(Form("%s", xsec_types_pretty.at(k_xsec_dataxsec).c_str() ));
    h_universe.at(k_up).at(k_xsec_dataxsec)->GetXaxis()->SetTitle("");
    h_universe.at(k_up).at(k_xsec_dataxsec)->GetXaxis()->SetLabelSize(0);



    h_universe.at(k_up).at(k_xsec_dataxsec)->Draw("hist");
    h_universe.at(k_dn).at(k_xsec_dataxsec)->Draw("hist,same");
    h_universe.at(k_up).at(k_xsec_mcxsec)->Draw("hist,same");
    h_universe.at(k_dn).at(k_xsec_mcxsec)->Draw("hist,same");


    // FIxed scaling for differential cross section
    if (vars.at(var) != "integrated"){
        // h_universe.at(k_up).at(k_xsec_dataxsec)->GetYaxis()->SetRangeUser(0, 0.5e-39);
    }

    TLegend *leg = new TLegend(0.5, 0.7, 0.85, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_universe.at(k_up).at(k_xsec_dataxsec), Form("%s +1 #sigma Data", label_pretty.c_str()), "l");
    leg->AddEntry(h_universe.at(k_up).at(k_xsec_mcxsec),   Form("%s +1 #sigma MC",   label_pretty.c_str()), "l");
    leg->AddEntry(h_universe.at(k_dn).at(k_xsec_dataxsec), Form("%s -1 #sigma Data", label_pretty.c_str()), "l");
    leg->AddEntry(h_universe.at(k_dn).at(k_xsec_mcxsec),   Form("%s -1 #sigma MC",   label_pretty.c_str()), "l");
    leg->Draw();

    bottomPad->cd();
    
    // Up ratio to CV
    TH1D* h_err_up = (TH1D *)h_universe.at(k_up).at(k_xsec_dataxsec)->Clone("h_ratio_up");
    h_err_up->Add(h_universe.at(k_up).at(k_xsec_mcxsec), -1);
    h_err_up->Divide(h_universe.at(k_up).at(k_xsec_dataxsec));
    h_err_up->SetLineWidth(2);
    h_err_up->SetLineColor(kGreen+2);
    h_err_up->Scale(100);
    
    // Down ratio to CV
    TH1D* h_err_dn = (TH1D *)h_universe.at(k_dn).at(k_xsec_dataxsec)->Clone("h_ratio_dn");
    h_err_dn->Add(h_universe.at(k_dn).at(k_xsec_mcxsec), -1);
    h_err_dn->Divide(h_universe.at(k_dn).at(k_xsec_dataxsec));
    h_err_dn->SetLineWidth(2);
    h_err_dn->SetLineColor(kRed+2);
    h_err_dn->Scale(100);

    TH1D* h_err = (TH1D *)cv_hist_vec.at(var).at(k_xsec_dataxsec)->Clone("h_ratio");
    h_err->Divide(cv_hist_vec.at(var).at(k_xsec_dataxsec));

    SetRatioOptions(h_err_up);
    h_err_up->GetYaxis()->SetNdivisions(4, 0, 0, kFALSE);
    h_err_up->GetYaxis()->SetRangeUser(-100, 100);
    h_err_up->GetYaxis()->SetTitle("\% change of Data to MC");
    h_err_up->GetXaxis()->SetTitle(var_labels_x.at(var).c_str());
    h_err_up->Draw("hist,same");
    h_err_dn->Draw("hist,same");
    h_err->Draw("hist,same");

    c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_data_mc_comparison.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(), _util.run_period, label.c_str(), vars.at(var).c_str() ));

    delete c;
}
// -----------------------------------------------------------------------------
void SystematicsHelper::CalcCovariance(std::string label, int var, std::vector<std::vector<TH1D*>> h_universe ){

    int n_bins = cv_hist_vec.at(var).at(0)->GetNbinsX();

    TH2D* cov  = new TH2D("h_cov",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the covariance matrix
    TH2D* cor  = new TH2D("h_cor",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the correlation matrix


    // Loop over universes
    std::cout << "Universes: " << h_universe.size() << std::endl;
    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over the rows
        for (int row = 1; row < cv_hist_vec.at(var).at(0)->GetNbinsX()+1; row++){
            
            double uni_row = h_universe.at(uni).at(k_xsec_mcxsec)->GetBinContent(row);
            double cv_row  = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(row);

            // Loop over the columns
            for (int col = 1; col < cv_hist_vec.at(var).at(0)->GetNbinsX()+1; col++){

                double uni_col = h_universe.at(uni).at(k_xsec_mcxsec)->GetBinContent(col);
                double cv_col  = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(col);
                
                double c = (uni_row - cv_row) * (uni_col - cv_col);

                if (uni != h_universe.size()-1)    cov->SetBinContent(row, col, cov->GetBinContent(row, col) + c ); // Fill with variance 
                else cov->SetBinContent(row, col, (cov->GetBinContent(row, col) + c) / h_universe.size());       // Fill with variance and divide by nuni
            
            }

        }
            
    }


    // ------------ Now calculate the correlation matrix ------------
    double cor_bini;
    // loop over rows
    for (int i=1; i<cv_hist_vec.at(var).at(0)->GetNbinsX()+1; i++) {
        double cii = cov->GetBinContent(i, i);

        // Loop over columns
        for (int j=1; j<cv_hist_vec.at(var).at(0)->GetNbinsX()+1; j++) {
            double cjj = cov->GetBinContent(j, j);
            double n = sqrt(cii * cjj);

            // Catch Zeros, set to arbitary 1.0
            if (n == 0) cor_bini = 0;
            else cor_bini = cov->GetBinContent(i, j) / n;

            cor->SetBinContent(i, j, cor_bini );
        }
    }

    TH2D *frac_cov = (TH2D*) cov->Clone("h_frac_cov");

    // ------------ Now calculate the fractional covariance matrix ------------
    double setbin;
    // loop over rows
    for (int i=1; i<cv_hist_vec.at(var).at(0)->GetNbinsX()+1; i++) {
        double cii = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(i);

        // Loop over columns
        for (int j=1; j<cv_hist_vec.at(var).at(0)->GetNbinsX()+1; j++) {
            double cjj = cv_hist_vec.at(var).at(k_xsec_mcxsec)->GetBinContent(j);
            double n = cii * cjj;

            // Catch Zeros, set to arbitary 0
            if (n == 0) setbin = 0;
            else setbin = frac_cov->GetBinContent(i, j) / n;

            frac_cov->SetBinContent(i, j, setbin );
        }
    }

    const Int_t NCont = 100;
    const Int_t NRGBs = 5;
    Double_t mainColour[NRGBs]   = { 1.00, 1.00, 1.00, 1.00, 1.00 };
    Double_t otherColour[NRGBs]   = { 0.99,0.80, 0.60, 0.40, 0.20 };
    Double_t stops[NRGBs] = { 0.00, 0.05, 0.1, 0.4, 1.00 };

    TColor::CreateGradientColorTable(NRGBs, stops, mainColour, otherColour, otherColour, NCont);
    gStyle->SetNumberContours(NCont);



    TCanvas *c = new TCanvas("c", "c", 500, 500);
    cov->GetXaxis()->CenterLabels(kTRUE);
    cov->GetYaxis()->CenterLabels(kTRUE);
    cov->SetTitle("Covariance Matrix");
    cov->Draw("colz");
    _util.IncreaseLabelSize(cov, c);
    c->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_cov.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str()));

    TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
    cor->GetXaxis()->CenterLabels(kTRUE);
    cor->GetYaxis()->CenterLabels(kTRUE);
    cor->SetTitle("Correlation Matrix");
    // cor->SetMaximum(1);
    // cor->SetMinimum(0);
    cor->Draw("colz, text00");
    _util.IncreaseLabelSize(cor, c2);
    cor->SetMarkerSize(1.3);
    gStyle->SetPaintTextFormat("0.3f");
    c2->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_cor.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str()));

    TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);
    frac_cov->GetXaxis()->CenterLabels(kTRUE);
    frac_cov->GetYaxis()->CenterLabels(kTRUE);
    frac_cov->Draw("colz, text00");
    frac_cov->SetTitle("Fractional Covariance Matrix");
    _util.IncreaseLabelSize(frac_cov, c3);
    frac_cov->SetMarkerSize(1.3);
    gStyle->SetPaintTextFormat("0.3f");
    c3->Print(Form("plots/run%s/Systematics/%s/%s/run%s_%s_%s_frac_cov.pdf", _util.run_period, label.c_str(), vars.at(var).c_str(),  _util.run_period, label.c_str(), vars.at(var).c_str()));

    delete cov;
    delete c;
    delete c2;
    delete cor;
    delete c3;
    delete frac_cov;

}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillSysVector(std::string variation, int var, int type, TH1D *h_up, TH1D *h_dn){

    // This is a Genie Unisim
    if (variation == "RPA"              ||
        variation == "CCMEC"            ||
        variation == "AxFFCCQE"         ||
        variation == "VecFFCCQE"        ||
        variation == "DecayAngMEC"      ||
        variation == "ThetaDelta2Npi"   ||
        variation == "ThetaDelta2NRad"  ||
        variation == "RPA_CCQE_Reduced" ||
        variation == "NormCCCOH"        ||
        variation == "NormNCCOH"){

        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double max_err = 0;
            // Get the max error in each bin, then add the square
            if (h_up->GetBinContent(bin+1) > h_dn->GetBinContent(bin+1)) max_err = h_up->GetBinContent(bin+1);
            else max_err = h_dn->GetBinContent(bin+1);
            v_genie_uni_total.at(var).at(type).at(bin) += max_err*max_err;
            
        }
    }
    // This is a beamline unisim
    else if (
        variation == "Horn_curr" ||
        variation == "Horn1_x" ||
        variation == "Horn1_y" ||
        variation == "Beam_spot" ||
        variation == "Horn2_x" ||
        variation == "Horn2_y" ||
        variation == "Horn_Water" ||
        variation == "Beam_shift_x" ||
        variation == "Beam_shift_y" ||
        variation == "Target_z" ||
        variation == "Horn1_refined_descr" ||
        variation == "Decay_pipe_Bfield" ||
        variation == "Old_Horn_Geometry"){

        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double max_err = 0;
            // Get the max error in each bin, then add the square
            if (h_up->GetBinContent(bin+1) > h_dn->GetBinContent(bin+1)) max_err = h_up->GetBinContent(bin+1);
            else max_err = h_dn->GetBinContent(bin+1);
            v_beamline_total.at(var).at(type).at(bin) += max_err*max_err;
            v_sys_total.at(var).at(type).at(bin) += max_err*max_err;
            
        }
    }
    else if (variation == "weightsGenie"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double max_err = 0;
            // Get the max error in each bin, then add the square
            if (h_up->GetBinContent(bin+1) > h_dn->GetBinContent(bin+1)) max_err = h_up->GetBinContent(bin+1);
            else max_err = h_dn->GetBinContent(bin+1);
            v_genie_multi_total.at(var).at(type).at(bin) += max_err*max_err;
            v_sys_total.at(var).at(type).at(bin) += max_err*max_err;
            
        }

    }
    else if (variation == "weightsReint"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double max_err = 0;
            // Get the max error in each bin, then add the square
            if (h_up->GetBinContent(bin+1) > h_dn->GetBinContent(bin+1)) max_err = h_up->GetBinContent(bin+1);
            else max_err = h_dn->GetBinContent(bin+1);
            v_reint_total.at(var).at(type).at(bin) += max_err*max_err;
            v_sys_total.at(var).at(type).at(bin) += max_err*max_err;
            
        }

    }
    else if (variation == "weightsPPFX"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double max_err = 0;
            // Get the max error in each bin, then add the square
            if (h_up->GetBinContent(bin+1) > h_dn->GetBinContent(bin+1)) max_err = h_up->GetBinContent(bin+1);
            else max_err = h_dn->GetBinContent(bin+1);
            v_hp_total.at(var).at(type).at(bin) += max_err*max_err;
            v_sys_total.at(var).at(type).at(bin) += max_err*max_err;
            
        }

    }
    else if (variation == "Dirt"){
        // Loop over histogram bins
        for (int bin = 0; bin < h_up->GetNbinsX(); bin++){
            
            double max_err = 0;
            // Get the max error in each bin, then add the square
            if (h_up->GetBinContent(bin+1) > h_dn->GetBinContent(bin+1)) max_err = h_up->GetBinContent(bin+1);
            else max_err = h_dn->GetBinContent(bin+1);
            v_dirt_total.at(var).at(type).at(bin) += max_err*max_err;
            v_sys_total.at(var).at(type).at(bin) += max_err*max_err;
            
        }

    }
    else {
        std::cout << "Unknown variation specified: " << variation << std::endl;
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillStatVector(){

    // Loop over the differential variables
    for (unsigned int var = 0; var < vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            
            // Get the uncertainty in each bin
            for (int bin = 0; bin < cv_hist_vec.at(var).at(type)->GetNbinsX(); bin++){
            
                double stat_err = 100 * cv_hist_vec.at(var).at(type)->GetBinError(bin+1) / cv_hist_vec.at(var).at(type)->GetBinContent(bin+1);
                v_stat_total.at(var).at(type).at(bin) += stat_err*stat_err;

            }
        
        }
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::FillPOTCountingVector(){

    // Loop over the differential variables
    for (unsigned int var = 0; var < vars.size(); var++ ){
        
        // Loop over the types
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            
            // Get the uncertainty in each bin
            for (int bin = 0; bin < cv_hist_vec.at(var).at(type)->GetNbinsX(); bin++){
            
                v_pot_total.at(var).at(type).at(bin) += 2.0*2.0;
                v_sys_total.at(var).at(type).at(bin) += 2.0*2.0;

            }
        
        }
    }

}
// -----------------------------------------------------------------------------
void SystematicsHelper::PrintUncertaintySummary(){

    // Loop over the variables
    for (unsigned int var = 0; var < vars.size(); var++ ){
        if (vars.at(var) == "true_el_E") continue; // Skip the true var which doesnt make much sense

        std::cout <<"----------------------------------------------" << std::endl;
        std::cout <<"Differential Variable: " << vars.at(var) <<"\n"<< std::endl;
        
        // Loop over the bins
        for (unsigned int bin = 0; bin < v_genie_uni_total.at(var).at(k_xsec_mcxsec).size(); bin++ ){
            
            std::cout << "Bin: " << bin+1 << " GENIE Unisim:       " <<std::sqrt(v_genie_uni_total.at(var).at(k_xsec_mcxsec)  .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " GENIE Multisim:     " <<std::sqrt(v_genie_multi_total.at(var).at(k_xsec_mcxsec).at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Beamline:           " <<std::sqrt(v_beamline_total.at(var).at(k_xsec_mcxsec)   .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Hadron Prod.:       " <<std::sqrt(v_hp_total.at(var).at(k_xsec_mcxsec)         .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Geant Rein.:        " <<std::sqrt(v_reint_total.at(var).at(k_xsec_mcxsec)      .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Dirt:               " <<std::sqrt(v_dirt_total.at(var).at(k_xsec_mcxsec)       .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " POT Counting:       " <<std::sqrt(v_pot_total.at(var).at(k_xsec_mcxsec)        .at(bin)) << " \%"<< std::endl;
            std::cout << std::endl;
            std::cout << "Bin: " << bin+1 << " Tot Data X-Sec Stat:                 " <<std::sqrt(v_stat_total.at(var).at(k_xsec_dataxsec)   .at(bin)) << " \%"<< std::endl;
            std::cout << "Bin: " << bin+1 << " Tot MC X-Sec Sys:                    " <<std::sqrt(v_sys_total.at(var).at(k_xsec_mcxsec)        .at(bin)) << " \%"<< std::endl;
            std::cout <<"\n--" << std::endl;       
        }
    }
    std::cout <<"----------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------
void SystematicsHelper::InitialiseUncertaintyVectors(){

    // differential variable, type, bin error
    v_genie_uni_total  .resize(vars.size());
    v_genie_multi_total.resize(vars.size());
    v_beamline_total   .resize(vars.size());
    v_hp_total         .resize(vars.size());
    v_reint_total      .resize(vars.size());
    v_sys_total        .resize(vars.size());
    v_stat_total       .resize(vars.size());
    v_dirt_total       .resize(vars.size());
    v_pot_total        .resize(vars.size());

    for (unsigned int var = 0; var < vars.size(); var++ ){
        v_genie_uni_total.at(var)  .resize(xsec_types.size());
        v_genie_multi_total.at(var).resize(xsec_types.size());
        v_beamline_total.at(var)   .resize(xsec_types.size());
        v_hp_total.at(var)         .resize(xsec_types.size());
        v_reint_total.at(var)      .resize(xsec_types.size());
        v_sys_total.at(var)        .resize(xsec_types.size());
        v_stat_total.at(var)       .resize(xsec_types.size());
        v_dirt_total.at(var)       .resize(xsec_types.size());
        v_pot_total .at(var)       .resize(xsec_types.size());
    }

    for (unsigned int var = 0; var < vars.size(); var++ ){
        for (unsigned int type = 0; type < xsec_types.size(); type++ ){
            v_genie_uni_total.at(var).at(type)  .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_genie_multi_total.at(var).at(type).resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_beamline_total.at(var).at(type)   .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_hp_total.at(var).at(type)         .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_reint_total.at(var).at(type)      .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_sys_total.at(var).at(type)        .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_stat_total.at(var).at(type)       .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_dirt_total.at(var).at(type)       .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
            v_pot_total .at(var).at(type)       .resize(cv_hist_vec.at(var).at(type)->GetNbinsX(), 0.0);
        }
    }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
