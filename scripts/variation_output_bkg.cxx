#include "variation_output_bkg.h"

//***************************************************************************
//***************************************************************************
int variation_output_bkg::GetLeadingShowerIndex(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
	int leading_index{0};
	int most_hits{0};
	
	// Loop over Particle Objects
	for (int j = 0; j < n_pfp; j++) {
		auto const pfp_obj = tpc_obj.GetParticle(j);
		const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
		
		if (pfp_pdg != 11) continue; // skip if not a shower

		const int n_pfp_hits = pfp_obj.NumPFPHits();
		
		// Compare for most hits
		if (n_pfp_hits > most_hits) {
			leading_index = j; 
			most_hits = n_pfp_hits; 
		}
	}

	return leading_index;
}
// Get the secondary shower index, if no secondary shower exists, returns -999
int variation_output_bkg::GetSecondaryShowerIndex(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
	int leading_index{0};
	int most_hits{0};
	int secondary_index{-999};
	
	// Loop over Particle Objects
	for (int j = 0; j < n_pfp; j++) {
		auto const pfp_obj = tpc_obj.GetParticle(j);
		const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
		
		if (pfp_pdg != 11) continue; // skip if not a shower

		const int n_pfp_hits = pfp_obj.NumPFPHits();
		
		// Compare for most hits
		if (n_pfp_hits > most_hits) {
			leading_index = j; 
			most_hits = n_pfp_hits; 
		}
	}
	// Loop over Particle Objects but this time omit leading index
	for (int j = 0; j < n_pfp; j++) {
		
		if(j == leading_index) continue;         //we assume leading shower == electron shower
		
		auto const pfp_obj = tpc_obj.GetParticle(j);
		const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
		
		if (pfp_pdg != 11) continue; // skip if not a shower if one doesnt exit returns -999

		const int n_pfp_hits = pfp_obj.NumPFPHits();
		
		// Compare for most hits
		if (n_pfp_hits > most_hits) {
			secondary_index = j; 
			most_hits = n_pfp_hits; 
		}
	}

	return secondary_index;
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::GetLongestTrackLength(const int n_pfp, int n_tpc_obj, xsecAna::TPCObjectContainer tpc_obj){
	double longest_track = 0;

	// Loop over Particle Objects
	for(int i = 0; i < n_pfp; i++) {

		auto const pfp_obj = tpc_obj.GetParticle(i);
		const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
		
		// If it  is track like
		if(pfp_pdg == 13) {
			const double trk_length = pfp_obj.pfpLength();
			
			if(trk_length > longest_track) {

				longest_track = trk_length;
			}
		}
	}

	return longest_track;
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::GetNumber_Track_Shower(const int n_pfp, int n_tpc_obj,
									 xsecAna::TPCObjectContainer tpc_obj, int &n_showers, int &n_tracks,
									 int &n_pfp_50Hits, int &n_tracks_50Hits,int &n_showers_50Hits) { 
	// Loop over Particle Objects
	for (int j = 0; j < n_pfp; j++) {
		auto const pfp_obj = tpc_obj.GetParticle(j);
		const int  pfp_pdg = pfp_obj.PFParticlePdgCode();
		
		// Don't count the neutrinos
		if (pfp_pdg != 12 || pfp_pdg != -12 || pfp_pdg != 14 || pfp_pdg != -14){
			// Cut with pfp with > 50 Hits
			if (pfp_obj.NumPFPHits() > 50){
				// For > 50 hits
				n_pfp_50Hits++; // Add one to n_pfp > 50 hits counter
				
				if (pfp_pdg == 11) 		n_showers_50Hits++; // Add to shower counter
				else if (pfp_pdg == 13) n_tracks_50Hits++;  // Add to track counter
				else std::cout << "Unknown pandora classification:\t" << pfp_pdg << std::endl;
					
			}

			// All PFP
			if (pfp_pdg == 11) 		n_showers++; // Add to shower counter
			else if (pfp_pdg == 13) n_tracks++;  // Add to track counter
			else return;
		}
		
	}
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::pfp_vtx_distance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z,
                                       double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z) {
	const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );

	return distance;
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::DrawTH1D(TH1D* h, double POT_Scaling){
	TCanvas* c = new TCanvas();
	c->cd();

	h->SetLineColor(kMagenta+3);
	h->SetLineWidth(2);
	h->SetLineStyle(1);

	h->Scale(POT_Scaling);
	h->Draw("his");

}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::DrawTH1D_SAME(TH1D* hist, std::string variation, TLegend* legend, std::string histname){
	
	// ----------------------
	//    Axis Specifiers
	// ----------------------
	if (histname == "h_total_hits"){
		hist->SetTitle(";PFP Total Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,4500);
	}
	else if (histname == "h_ldg_shwr_hits") {
		hist->SetTitle("; Leading Shower Hits (All Planes);Entries");
		// hist->GetYaxis()->SetRangeUser(0,750);
	}
	else if (histname == "h_ldg_shwr_hits_WPlane") {
		hist->SetTitle(";Leading Shower Collection Plane Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,5000);
	}
	else if (histname == "h_ldg_shwr_Open_Angle"){
		hist->SetTitle(";Leading Shower Opening Angle [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,175);
	}
	else if (histname == "h_ldg_shwr_dEdx_WPlane"){
		hist->SetTitle(";Leading Shower dEdx Collection Plane [MeV/cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_HitPerLen"){
		hist->SetTitle(";Leading Shower Hits / Length [ cm^{-1} ];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_Phi"){
		hist->SetTitle(";Leading Shower #phi [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_Theta"){
		hist->SetTitle(";Leading Shower #theta [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_ldg_shwr_CTheta"){
		hist->SetTitle(";Leading Shower cos(#theta);Entries");
		// hist->GetYaxis()->SetRangeUser(0,7000);
	}
	else if (histname == "h_long_Track_ldg_shwr"){
		hist->SetTitle(";Longest Track Length / Leading Shower Length;Entries");
		// hist->GetYaxis()->SetRangeUser(0,22000);
	}
	else if (histname == "h_tpc_obj_vtx_x"){
		hist->SetTitle(";TPC Object Vertex X [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,2500);
	}
	else if (histname == "h_tpc_obj_vtx_y"){
		hist->SetTitle(";TPC Object Vertex Y [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,2200);
	}
	else if (histname == "h_tpc_obj_vtx_z"){
		hist->SetTitle(";TPC Object Vertex Z [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_pfp"){
		hist->SetTitle(";Number of PFP;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_pfp_50Hits"){
		hist->SetTitle("; Number of PFP > 50 Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_tracks"){
		hist->SetTitle(";Number of Tracks;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_tracks_50Hits"){
		hist->SetTitle("; Number of Tracks > 50 Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_showers"){
		hist->SetTitle(";Number of Showers;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_n_showers_50Hits"){
		hist->SetTitle("; Number of Showers > 50 Hits;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_track_phi"){
		hist->SetTitle("; Track #phi [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,6000);
	}
	else if (histname == "h_shower_phi"){
		hist->SetTitle("; Shower #phi [degrees];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}			
	else if (histname == "h_largest_flash_y"){
		hist->SetTitle(";Largest Flash Y [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_largest_flash_z"){
		hist->SetTitle("; Largest Flash Z [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1000);
	}
	else if (histname == "h_largest_flash_time"){
		hist->SetTitle("; Largest Flash Time [#mus];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_largest_flash_pe"){
		hist->SetTitle("; Largest Flash PE;Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_Flash_TPCObj_Dist"){
		hist->SetTitle("; 2D Distance from Largest Flash to #nu Vertex [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_shower_Nu_vtx_Dist"){
		hist->SetTitle("; 3D Distance of Shower to #nu Vertex [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,1200);
	}
	else if (histname == "h_track_Nu_vtx_Dist"){
		hist->SetTitle("; 3D Distance of Track to #nu Vertex [cm];Entries");
		// hist->GetYaxis()->SetRangeUser(0,18000);
	}
	else return;

	// ----------------------
	//    Draw Specifiers
	// ----------------------
	if (variation == "BNBCV" || variation == "NuMICV"){
		hist->SetLineColor(kBlack);
		hist->SetLineWidth(2);
		hist->SetLineStyle(1);
		legend->AddEntry(hist, "CV", "l");
		hist->Draw("hist,same");
	} 
	else if  (variation == "BNBwithDIC"         || variation == "NuMIwithDIC"){
		hist->SetLineColor(kMagenta+2);
		hist->SetLineWidth(2);
		hist->SetLineStyle(1);
		legend->AddEntry(hist, "DIC", "l");
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBEnhancedTPCVis"  || variation == "NuMIEnhancedTPCVis" ){ 
		hist->SetLineColor(30);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Enhanced TPC Vis.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBaltDeadChannels" || variation == "NuMIaltDeadChannels"){ 
		hist->SetLineColor(38);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Alt. Dead Chan.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBdeadSaturatedChannels" || variation == "NuMIdeadSaturatedChannels"){
		hist->SetLineColor(28);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Dead Sat. Chan.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
		
	}
	else if  (variation == "BNBstretchResp"  || variation == "NuMIstretchResp" ){
		hist->SetLineColor(36);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Stretch Resp.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBsqueezeResp"  || variation == "NUMIsqueezeResp"){
		hist->SetLineColor(1001);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Squeeze Resp.", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBupPEnoise"    || variation == "NuMIupPEnoise"){
		hist->SetLineColor(kBlue+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "PE Noise Up", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBnoiseAmpDown" || variation == "NUMInoiseAmpDown"){
		hist->SetLineColor(42);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Noise Amp. Down", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBdownPEnoise"  || variation == "NUMIdownPEnoise"){
		hist->SetLineColor(50);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "PE Noise Down", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBnoiseAmpUp"   || variation == "NUMInoiseAmpUp"){
		hist->SetLineColor(kOrange+10);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Noise Amp. Up", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBDTdown"       || variation == "NUMIDTdown"){
		hist->SetLineColor(kOrange+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "DT Down", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBDTup"         || variation == "NUMIDTup"){
		hist->SetLineColor(kMagenta-10);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "DT Up", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBDLup"         || variation == "NUMIDLup"){
		hist->SetLineColor(kMagenta);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "DL Up", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBDLdown"       || variation == "NUMIDLdown"){
		hist->SetLineColor(kTeal+6);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "DL Down", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBdataSCE"      || variation == "NUMIdataSCE"){
		hist->SetLineColor(kYellow);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "SCE", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
	else if  (variation == "BNBLArG4BugFix"  || variation == "NUMILArG4BugFix"){
		hist->SetLineColor(kSpring-7);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "LArG4BugFix", "l");
		hist->SetLineStyle(1);
		hist->Draw("hist,same");
	}
    else if  (variation == "BNBBirksRecomb"  || variation == "NUMIBirksRecomb"){
        hist->SetLineColor(kRed+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Birks Recomb.","l");
        hist->SetLineStyle(1);
        hist->Draw("hist,same");
    }
	else {
		std::cout << "Error! Could not match varations:\t" << variation << std::endl;
		return;
	}
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::GetPOT(const char * _file1){
	double POT{0};
	std::string line;

	std::string filename;
	std::string temp_filename = _file1; // cast to string
	filename = "File: " + temp_filename;

	std::cout << filename << std::endl;

	std::ifstream myfile ("POT_List_Varaitions.txt");
	int i_POT{0};
	
	if (myfile.is_open()) {

		// Loop over lines in file
		while ( getline (myfile,line) ) {

			if (i_POT == 1){
				line.erase(0, 21);
				POT = std::stod(line); // Convert string to double
				std::cout << "POT in File:\t"<< POT << std::endl;
				return POT;
			}

			// Found the correct variation file 
			if (line == filename){
				std::cout << "Found match for POT file"<< std::endl;
				i_POT++;
			}
			
		}

		myfile.close();
	}

	else std::cout << "Unable to open file" << std::endl; 
	std::cout << "Could not find a match for POT file"<< std::endl;
	return POT;
}
//***************************************************************************
//***************************************************************************
// Function that grabs the reweighted histogram names for plotting
std::vector<std::string> variation_output_bkg::GrabDirs(TFile* f_var_out) {
	std::vector<std::string> variations;

	f_var_out->cd();
	
	TKey *key;
	TIter nextkey(gDirectory->GetListOfKeys());

	TString same_plots = "SAME_Plots";

	std::cout << "\n=================================================" << std::endl;	
	std::cout << "Getting variation modes:" << std::endl;	
  	while ( ( key =  (TKey*)nextkey()) ) { // Extra brackets to omit a warning 
		if (key->IsFolder()) {
			std::cout << key->GetName() << std::endl; // Print the variations
			if (key->GetName() == same_plots ) continue; // Skip this
			variations.push_back(key->GetName());
		}
	}
	std::cout << "=================================================\n" << std::endl;

	return (variations);
}
//***************************************************************************
//***************************************************************************
void variation_output_bkg::PlotVariatons(TFile* f_var_out){
	f_var_out->cd();
	
	// Grab the variation folders in the file
	std::vector<std::string> variations = variation_output_bkg::GrabDirs(f_var_out); 

	TH1D* hist;
	
	std::vector<std::string> histnames = {"h_total_hits","h_ldg_shwr_hits", "h_ldg_shwr_hits_WPlane",
										 "h_ldg_shwr_Open_Angle", "h_ldg_shwr_dEdx_WPlane", "h_ldg_shwr_HitPerLen",
										 "h_ldg_shwr_Phi", "h_ldg_shwr_Theta","h_ldg_shwr_CTheta",
										  "h_long_Track_ldg_shwr", "h_tpc_obj_vtx_x", "h_tpc_obj_vtx_y", "h_tpc_obj_vtx_z",
										  "h_n_pfp", "h_n_pfp_50Hits", "h_n_tracks", "h_n_tracks_50Hits", "h_n_showers",
										  "h_n_showers_50Hits", "h_track_phi", "h_shower_phi", "h_largest_flash_y", "h_largest_flash_z",
										  "h_largest_flash_time", "h_largest_flash_pe", "h_Flash_TPCObj_Dist",
										  "h_shower_Nu_vtx_Dist", "h_track_Nu_vtx_Dist" };

	// Loop over the histograms
	for (int j=0; j < histnames.size(); j++){
		
		// Canvas + Legend
		TCanvas* c = new TCanvas();
		TLegend* legend;
		if (histnames[j] == "h_ldg_shwr_CTheta") legend = new TLegend(0.15, 0.55, 0.37, 0.85); // Reposition
		else if (histnames[j] == "h_tpc_obj_vtx_x" || histnames[j] == "h_tpc_obj_vtx_y" || histnames[j] == "h_tpc_obj_vtx_z" )
				legend = new TLegend(0.15, 0.15, 0.37, 0.45); // Reposition
		else if (histnames[j] == "h_largest_flash_z" ) legend = new TLegend(0.35, 0.59, 0.57, 0.89);
		else legend = new TLegend(0.72, 0.59, 0.94, 0.89);

		legend->SetBorderSize(0);
		legend->SetFillStyle(0);

		// Loop over variation directories
		for (int i=0; i < variations.size(); i++){
			char name[500];
			snprintf(name, 500, "%s/%s", variations[i].c_str(),histnames[j].c_str() );

			hist = (TH1D*)f_var_out->Get(name);

			if (hist == NULL ) std::cout << "ERROR: Can't get Histogram!" << std::endl;

			DrawTH1D_SAME(hist, variations[i], legend, histnames[j]);

		}

		// Print the Canvas
		char Canvas_name[500];
		snprintf(Canvas_name, 500, "plots/%s.pdf",histnames[j].c_str() ); 
		legend->Draw();
		c->Print(Canvas_name);

	}

	// Close the file
	f_var_out->Close();

}
//***************************************************************************
//************************** Flash Functions ********************************
std::vector<std::vector<double>> variation_output_bkg::GetLargestFlashVector(TFile* f, double flash_start_time, double flash_end_time ){

	// Get Optical Information from file
	TTree* optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	std::vector<std::vector<double>> largest_flash_v_v;

	// ----------------------
	//    Optical Info
	// ----------------------
	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;

	optical_tree->SetBranchAddress("Run",              &fRun           ); // Run number of Flash
	optical_tree->SetBranchAddress("Event",            &fEvent         ); // Event number of Flash
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE     ); // PE of Flash
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime   ); // Time of flash
	optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY );
	optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ );
	optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY); // Flash Y center
	optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ); // Flash Z Center

	// Num events in the optical tree
	const int optical_entries = optical_tree->GetEntries();
	std::cout << "Total Optical Entries: " << optical_entries << std::endl;

	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;

	// Contains the entry number for a given OpFlash per event
	std::vector<int>					optical_list_pe;
	std::vector<std::vector<int> >		optical_list_pe_v;
	
	std::vector<double>					optical_list_flash_center_y; 
	std::vector<std::vector<double> >	optical_list_flash_center_y_v;
	
	std::vector<double>					optical_list_flash_center_z; 
	std::vector<std::vector<double> >	optical_list_flash_center_z_v;
	
	std::vector<double>					optical_list_flash_time;
	std::vector<std::vector<double> >	optical_list_flash_time_v;
	
	// Loop over the optical entries to get the largest flash vector
	
	// ----------------------
	// Resize the optical enties to be the same sizd as number of Events (TPC Obj)
	// ----------------------
	
	for(int i = 0; i < optical_entries; i++) {
		
		// Get the Optical entry
		optical_tree->GetEntry(i);

		current_run		= fRun;
		current_event 	= fEvent;

		// New event
		if(current_event != last_event) {
			optical_list_pe.clear();
			optical_list_flash_center_y.clear();
			optical_list_flash_center_z.clear();
			optical_list_flash_time.clear();

			optical_list_pe.push_back(fOpFlashPE);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);

		}
		// Same event
		if(current_event == last_event && current_run == last_run) {
			optical_list_pe_v.pop_back();
			optical_list_flash_center_y_v.pop_back();
			optical_list_flash_center_z_v.pop_back();
			optical_list_flash_time_v.pop_back();

			optical_list_pe.push_back(fOpFlashPE);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);

		}

		last_event = current_event;
		last_run   = current_run;

		optical_list_pe_v.push_back(optical_list_pe);
		optical_list_flash_center_y_v.push_back(optical_list_flash_center_y);
		optical_list_flash_center_z_v.push_back(optical_list_flash_center_z);
		optical_list_flash_time_v.push_back(optical_list_flash_time);

	}
	
	std::cout << "Resized Optical List Size: " << optical_list_pe_v.size() << std::endl;
	
	// Largest Flash Vector
	std::vector<double> largest_flash_v;
	largest_flash_v.resize(4, 0);
	
	// ----------------------
	//      Event loop
	// ----------------------
	for(int i = 0; i < optical_list_pe_v.size(); i++) {
		
		bool in_time 				= false;
		bool got_in_time 			= false;
		bool sufficient_flash 		= false;
		bool got_sufficient_flash 	= false;
		
		double largest_flash = 0.;
		double largest_center_y = 0;
		double largest_center_z = 0;
		double largest_flash_time = 0;
		
		// Cut Variables defined in main.h
		int flash_pe_threshold = 50;

		// Loop through all flashes in event and find largest
		for(int j = 0; j < optical_list_pe_v.at(i).size(); j++) {
			
			auto const opt_time         = optical_list_flash_time_v.at(i).at(j) + 0.94; // shift due to MC and offset
			auto const opt_pe           = optical_list_pe_v.at(i).at(j);
			const double opt_center_y   = optical_list_flash_center_y_v.at(i).at(j);
			const double opt_center_z   = optical_list_flash_center_z_v.at(i).at(j);
			const double opt_flash_time = optical_list_flash_time_v.at(i).at(j) + 0.94; // shift due to MC and offset
			
			// See if flash was in time
			in_time = flash_in_time(opt_time, flash_time_start, flash_time_end); 
			if(in_time == true) got_in_time = true; 
			
			// See if flash meets the threshold requirements
			sufficient_flash = flash_pe(opt_pe, flash_pe_threshold); 
			if(sufficient_flash == true) got_sufficient_flash = true; 
			
			//Flash is both in time and over PE threshold
			if(in_time == true && sufficient_flash == true) {
				
				// Find the largest flash in this event
				if(opt_pe > largest_flash) {
					largest_flash      = opt_pe;
					largest_center_y   = opt_center_y;
					largest_center_z   = opt_center_z;
					largest_flash_time = opt_flash_time;
				}
			}
		}
		
		largest_flash_v.at(0) = largest_center_y;
		largest_flash_v.at(1) = largest_center_z;
		largest_flash_v.at(2) = largest_flash_time;
		largest_flash_v.at(3) = largest_flash;
		
		largest_flash_v_v.push_back(largest_flash_v);
		
	}

	return largest_flash_v_v;
}
//***************************************************************************
//***************************************************************************
double variation_output_bkg::Flash_TPCObj_vtx_Dist(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z) {
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
	return distance;
}
//***************************************************************************
//***************************************************************************
std::pair<std::string, int> variation_output_bkg::TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool true_in_tpc, bool has_pi0) {
	int part_nue_cc     = 0;
	int part_nue_bar_cc = 0;
	int part_cosmic     = 0;
	int part_nc         = 0;
	int part_nc_pi0     = 0;
	int part_numu_cc    = 0;
	int part_unmatched  = 0;

	const int tpc_obj_mode = tpc_obj.Mode();
	const int n_pfp = tpc_obj.NumPFParticles();
	const int n_pfp_showers = tpc_obj.NPfpShowers();
	int most_hits = 0;
	int leading_index = -1;
	int leading_pdg = 0;
	int leading_mc_parent_pdg = 0;
	std::string leading_origin = "kNothing";

	for(int j = 0; j < n_pfp; j++)
	{
		auto const part = tpc_obj.GetParticle(j);
		const int n_pfp_hits = part.NumPFPHits();
		const int mc_parent_pdg = part.MCParentPdg();
		const int pfp_pdg = part.PFParticlePdgCode();
		if(pfp_pdg == 11)
		{
			if(n_pfp_hits > most_hits)
			{
				leading_index = j;
				most_hits = n_pfp_hits;
			}
		}
		//if(n_pfp_showers)
		if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && mc_parent_pdg == 12)  { part_nue_cc++; }
		if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && mc_parent_pdg == -12) { part_nue_bar_cc++; }
		if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (mc_parent_pdg == 14 || mc_parent_pdg == -14)) { part_numu_cc++; }
		if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino")
		{
			if(has_pi0 == true)  {part_nc_pi0++; }
			if(has_pi0 == false) {part_nc++; }
		}
		if(part.Origin() == "kCosmicRay") { part_cosmic++;    }
		if(part.Origin() == "kUnknown"  ) { part_unmatched++; }
	}
	//some tpc objects actually have 0 hits - crazy!
	if(tpc_obj.NumPFPHits() == 0) {return std::make_pair("bad_reco", 0); }

	//currently, any tpc objects which only have a track end up with a leading_index of -1
	//this index will likely cause code to crash if called before the signal definition cuts

	//also some rare cases where nu_pfp = nue, and shower hits = 0 with track hits > 0 - how does this happen? (NC event?)

	//now to catagorise the tpco
	if(part_cosmic > 0)
	{
		if(part_nue_cc  > 0 || part_nue_bar_cc > 0)  { return std::make_pair("nue_cc_mixed",  leading_index); }
		if(part_numu_cc > 0 )                        { return std::make_pair("numu_cc_mixed", leading_index); }
		if(part_nc  > 0 || part_nc_pi0 > 0)          { return std::make_pair("other_mixed",   leading_index); }
		return std::make_pair("cosmic", leading_index);
	}
	//this uses the true neutrino vertex for this specific event
	//not the true vtx per tpc object - maybe this can be fixed in the future...
	//but using the true nu vtx only matters for the pure signal events,
	//where the neutrino vertex IS the true tpc object vertex
	if(part_cosmic == 0)
	{
		if(part_nue_cc      > 0 && true_in_tpc == false) { return std::make_pair("nue_cc_out_fv", leading_index);   }
		if(part_nue_bar_cc  > 0 && true_in_tpc == false) { return std::make_pair("nue_cc_out_fv", leading_index);   }

		if(part_nue_cc    > 0 && tpc_obj_mode == 0   ) { return std::make_pair("nue_cc_qe",     leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 1   ) { return std::make_pair("nue_cc_res",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 2   ) { return std::make_pair("nue_cc_dis",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 3   ) { return std::make_pair("nue_cc_coh",    leading_index);   }
		if(part_nue_cc    > 0 && tpc_obj_mode == 10  ) { return std::make_pair("nue_cc_mec",    leading_index);   }

		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 0   ) { return std::make_pair("nue_bar_cc_qe",     leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 1   ) { return std::make_pair("nue_bar_cc_res",    leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 2   ) { return std::make_pair("nue_bar_cc_dis",    leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 3   ) { return std::make_pair("nue_bar_cc_coh",    leading_index);   }
		if(part_nue_bar_cc    > 0 && tpc_obj_mode == 10  ) { return std::make_pair("nue_bar_cc_mec",    leading_index);   }

		if(part_numu_cc     > 0 && tpc_obj_mode == 0   ) { return std::make_pair("numu_cc_qe",    leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 1   ) { return std::make_pair("numu_cc_res",   leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 2   ) { return std::make_pair("numu_cc_dis",   leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 3   ) { return std::make_pair("numu_cc_coh",   leading_index);   }
		if(part_numu_cc     > 0 && tpc_obj_mode == 10  ) { return std::make_pair("numu_cc_mec",   leading_index);   }
		if(part_nc          > 0                        ) { return std::make_pair("nc",            leading_index);   }
		if(part_nc_pi0      > 0                        ) { return std::make_pair("nc_pi0",        leading_index);   }
		if(part_unmatched   > 0                        ) { return std::make_pair("unmatched",     leading_index);   }
	}
	//this never happens :)
	std::cout << "HELP HELP HELP END OF TPCO CLASSIFIER AND NO CLASSIFICATION!" << std::endl;
	//return the string for the tpco id
	return std::make_pair("non_match", 0);
}
//*********************** Selection Cuts ************************************
//***************************************************************************
//IN FV
bool variation_output_bkg::in_fv(double x, double y, double z, std::vector<double> fv_boundary_v) {
	const double det_x1 = 0;
	const double det_x2 = 256.35;
	const double det_y1 = -116.5;
	const double det_y2 = 116.5;
	const double det_z1 = 0;
	const double det_z2 = 1036.8;

	const double x1 = fv_boundary_v.at(0);
	const double x2 = fv_boundary_v.at(1);
	const double y1 = fv_boundary_v.at(2);
	const double y2 = fv_boundary_v.at(3);
	const double z1 = fv_boundary_v.at(4);
	const double z2 = fv_boundary_v.at(5);

	if(x <= det_x1 + x1 || x >= det_x2 - x2) {return false; }
	if(y <= det_y1 + y1 || y >= det_y2 - y2) {return false; }
	if(z <= det_z1 + z1 || z >= det_z2 - z2) {return false; }
	return true;
}
//***************************************************************************
//***************************************************************************
// Flash in time
bool variation_output_bkg::flash_in_time(double flash_time, double flash_start, double flash_end) {
	if(flash_time >= flash_start && flash_time <= flash_end) return true; // Pass in time
	else return false; // Fail in time
}
//***************************************************************************
//***************************************************************************
// Flash PE
bool variation_output_bkg::flash_pe(int flash_pe, int flash_pe_threshold) {
	if (flash_pe >= flash_pe_threshold) return true; // Pass PE Thresh
	else return false; // Fail PE Thresh
}
//***************************************************************************
//***************************************************************************
// Get the vector of flashes with true/false on whether it passed the selection or not
void variation_output_bkg::FlashinTime_FlashPE(TFile* f, double flash_start_time, double flash_end_time, std::vector<bool> &flash_cuts_pass_vec ){
	// Get Optical Information from file
	TTree* optical_tree = (TTree*)f->Get("AnalyzeTPCO/optical_tree");

	std::vector<std::vector<double>> largest_flash_v_v;

	// ----------------------
	//    Optical Info
	// ----------------------
	int fRun = 0;
	int fEvent = 0;
	int fOpFlashPE = 0;
	double fOpFlashTime = 0;
	double fOpFlashWidthY = 0;
	double fOpFlashWidthZ = 0;
	double fOpFlashCenterY = 0;
	double fOpFlashCenterZ = 0;

	optical_tree->SetBranchAddress("Run",              &fRun           ); // Run number of Flash
	optical_tree->SetBranchAddress("Event",            &fEvent         ); // Event number of Flash
	optical_tree->SetBranchAddress("OpFlashPE",        &fOpFlashPE     ); // PE of Flash
	optical_tree->SetBranchAddress("OpFlashTime",      &fOpFlashTime   ); // Time of flash
	optical_tree->SetBranchAddress("OpFlashWidhtY",    &fOpFlashWidthY );
	optical_tree->SetBranchAddress("OpFlashWidthZ",    &fOpFlashWidthZ );
	optical_tree->SetBranchAddress("OpFlashCenterY",   &fOpFlashCenterY); // Flash Y center
	optical_tree->SetBranchAddress("OpFlashCenterZ",   &fOpFlashCenterZ); // Flash Z Center

	// Num events in the optical tree
	const int optical_entries = optical_tree->GetEntries();
	std::cout << "Total Optical Entries: " << optical_entries << std::endl;

	int current_event = 0;
	int current_run = 0;
	int last_event = 0;
	int last_run = 0;

	// Contains the entry number for a given OpFlash per event
	std::vector<int>					optical_list_pe;
	std::vector<std::vector<int> >		optical_list_pe_v;
	
	std::vector<double>					optical_list_flash_center_y; 
	std::vector<std::vector<double> >	optical_list_flash_center_y_v;
	
	std::vector<double>					optical_list_flash_center_z; 
	std::vector<std::vector<double> >	optical_list_flash_center_z_v;
	
	std::vector<double>					optical_list_flash_time;
	std::vector<std::vector<double> >	optical_list_flash_time_v;
	
	// Loop over the optical entries to get the largest flash vector
	
	// ----------------------
	// Resize the optical enties to be the same sizd as number of Events (TPC Obj)
	// ----------------------
	
	for(int i = 0; i < optical_entries; i++) {
		
		// Get the Optical entry
		optical_tree->GetEntry(i);

		current_run		= fRun;
		current_event 	= fEvent;

		// New event
		if(current_event != last_event) {
			optical_list_pe.clear();
			optical_list_flash_center_y.clear();
			optical_list_flash_center_z.clear();
			optical_list_flash_time.clear();

			optical_list_pe.push_back(fOpFlashPE);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);

		}
		// Same event
		if(current_event == last_event && current_run == last_run) {
			optical_list_pe_v.pop_back();
			optical_list_flash_center_y_v.pop_back();
			optical_list_flash_center_z_v.pop_back();
			optical_list_flash_time_v.pop_back();

			optical_list_pe.push_back(fOpFlashPE);
			optical_list_flash_center_y.push_back(fOpFlashCenterY);
			optical_list_flash_center_z.push_back(fOpFlashCenterZ);
			optical_list_flash_time.push_back(fOpFlashTime);

		}

		last_event = current_event;
		last_run   = current_run;

		optical_list_pe_v.push_back(optical_list_pe);
		optical_list_flash_center_y_v.push_back(optical_list_flash_center_y);
		optical_list_flash_center_z_v.push_back(optical_list_flash_center_z);
		optical_list_flash_time_v.push_back(optical_list_flash_time);

	}
	
	std::cout << "Resized Optical List Size: " << optical_list_pe_v.size() << std::endl;
	
	flash_cuts_pass_vec.resize(optical_list_pe_v.size());

	// Loop over the optical list
	for(int i = 0; i < optical_list_pe_v.size(); i++) {
		
		bool in_time = false;
		bool sufficient_flash = false;
		
		auto const opt_time_v = optical_list_flash_time_v.at(i);
		auto const opt_pe_v   = optical_list_pe_v.at(i);
		
		// Loop over the optical list vec
		for (int j = 0; j < optical_list_pe_v.at(i).size(); j++) {
			
			auto const opt_time = opt_time_v.at(j) +1.0;
			auto const opt_pe = opt_pe_v.at(j);
			
			in_time          = flash_in_time(opt_time, flash_start_time, flash_end_time);
			sufficient_flash = flash_pe(opt_pe, flash_pe_threshold);
			
			// Flash is both in time and over PE threshold
			if(in_time == true && sufficient_flash == true){
				flash_cuts_pass_vec.at(i) = true;
				break; // once pased we are done, so dont loop any more otherwise we may overwrite this

			}
		}

		if (in_time == false && sufficient_flash == false) {
			flash_cuts_pass_vec.at(i) = false;
		}

	}

}
//***************************************************************************
//***************************************************************************
// Check for valid reconstructed nue
bool variation_output_bkg::HasNue(xsecAna::TPCObjectContainer tpc_obj, const int n_pfp ) {

	bool has_nue = false;
	bool has_valid_shower = false;

	for(int j = 0; j < n_pfp; j++) {
			auto const part     = tpc_obj.GetParticle(j);
			const int  pfp_pdg  = part.PFParticlePdgCode();
			const int  pfp_hits = part.NumPFPHits();
			
			if(pfp_pdg == 11 && pfp_hits > 0) has_valid_shower = true; 
			
			if(pfp_pdg == 12) has_nue = true; 
	}

	
	if(has_nue == true && has_valid_shower == true)
		return true; 
	else 
		return false;
	
}
//***************************************************************************
//***************************************************************************
// Flash reco Vertex Distance 
bool variation_output_bkg::opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance) {
	const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
	
	if(distance <= tolerance) return true;
	return false;
}
void variation_output_bkg::flashRecoVtxDist(std::vector< double > largest_flash_v, double tolerance, const double tpc_vtx_x, const double tpc_vtx_y, const double tpc_vtx_z, std::vector<bool> &flash_cuts_pass_vec) {
	
	bool is_close;
	for (unsigned int i=0; i < flash_cuts_pass_vec.size(); i++){
		
		//flash is upstream
		if(tpc_vtx_z < largest_flash_v.at(1)) 
			is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tolerance);
		
		//flash is downstream
		if(tpc_vtx_z >= largest_flash_v.at(1)) 
			is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), (tolerance - 20));
		
		
		if(is_close == true )
			flash_cuts_pass_vec.at(i) = true;
		
		if(is_close == false) 
			flash_cuts_pass_vec.at(i) = false;

	}

}
//***************************************************************************
//***************************************************************************
// Vertex to shower distance
bool variation_output_bkg::VtxNuDistance(double tpc_vtx_x, double tpc_vtx_y, double tpc_vtx_z, double pfp_vtx_x, double pfp_vtx_y, double pfp_vtx_z, double tolerance){
	
	const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );
	
	if(distance <= tolerance)
		return true; 
	else 
		return false;

}
//***************************************************************************
//***************************************************************************
// Hit thresholds all planes
bool variation_output_bkg::HitThreshold(const int n_pfp_hits, double threshold){

	if(n_pfp_hits >= threshold)
		return true;
	else
		return false; 

}
//***************************************************************************
//***************************************************************************
// Leading shower open angle
bool variation_output_bkg::OpenAngleCut(const double leading_open_angle, const std::vector<double> tolerance_open_angle){

	if(leading_open_angle <= tolerance_open_angle.at(1) && leading_open_angle >= tolerance_open_angle.at(0))
		return true;
	else 
		return false;

}
//***************************************************************************
//***************************************************************************
// leading shower dEdx cut
bool variation_output_bkg::dEdxCut( const double leading_dedx, const double tolerance_dedx_min, const double tolerance_dedx_max){

	if(leading_dedx <= tolerance_dedx_max && leading_dedx >= tolerance_dedx_min)
		return true;
	else 
		return false;

}
//***************************************************************************
//***************************************************************************
bool variation_output_bkg::HitLengthRatioCut(const double pfp_hits_length_tolerance, const double pfp_hits, const double pfp_length){
	const double pfp_hits_length_ratio = (pfp_hits / pfp_length);

	if (pfp_hits_length_ratio > pfp_hits_length_tolerance)
		return true;
	else
		return false;
	
}
//***************************************************************************
//***************************************************************************
bool variation_output_bkg::LongestTrackLeadingShowerCut(const double ratio_tolerance, double longest_track_leading_shower_ratio){

	if(longest_track_leading_shower_ratio < ratio_tolerance)
		return true;
	else
		return false;

}
//***************************************************************************
//***************************************************************************
bool variation_output_bkg::IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v) {
	
	if(in_fv(track_start.at(0), track_start.at(1), track_start.at(2), fv_boundary_v) == true
	   && in_fv(track_end.at(0), track_end.at(1), track_end.at(2), fv_boundary_v) == true) {
		return true;
	}
	else 
		return false;
}
bool variation_output_bkg::ContainedTracksCut(std::vector<double> fv_boundary_v, std::vector<double> pfp_start_vtx, std::vector<double> pfp_end_vtx){

	const bool is_contained = IsContained(pfp_start_vtx, pfp_end_vtx, fv_boundary_v);

	//if not contained
	if(is_contained == true) 
		return true;
	else
		return false;
	
}
//***************************************************************************
//***************************************************************************
// ----------------------
//         Main
// ----------------------
void variation_output_bkg::run_var(const char * _file1, TString mode, const std::vector<double> _config) {
	// std::cout << "=================================================\n" << std::endl;
	// std::cout << "Warning, there are hardcoded values" << 
	// 	"in this script, grep for  \"HARDCODED\" for places\n" << std::endl;
	// std::cout << "=================================================\n" << std::endl;
	
	// create plots folder if it does not exist
	system("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi");

	gStyle->SetOptStat(0); // say no to stats box

	// Choose while file 
	if (mode == "bnb")
		f_var_out = new TFile("plots/variation_out_bnb_bkg.root","UPDATE");
	else if (mode == "numi")
		f_var_out = new TFile("plots/variation_out_numi_bkg.root","UPDATE");
	else f_var_out = new TFile("plots/variation_out_bnb_bkg.root","UPDATE");

	//*************************** Configure the cut parameters *************************************
	std::cout << "\n --- Configuring Parameters --- \n" << std::endl;
	_x1 = _config[0];
	_x2 = _config[1];
	_y1 = _config[2];
	_y2 = _config[3];
	_z1 = _config[4];
	_z2 = _config[5];
	flash_pe_threshold = _config[6];
	
	if (mode == "bnb"){
		std::cout << "Using BNB Params" << std::endl;
		flash_time_start = 3.57; // Manually Override for BNB [us]
		flash_time_end   = 5.25; // Manually Override for BNB [us]
	}
	else if (mode == "numi") {
		std::cout << "Using NuMI Params" << std::endl;
		flash_time_start = _config[7]; // Use default numi config
		flash_time_end   = _config[8]; // Use default numi config
	}
	else if (mode == "same"){
		std::cout << "Forget about the configs, we are going to make some plots!" << std::endl;
	}
	else {
		std::cout << "Using default numi flash config since no match could be made"<< std::endl;
		flash_time_start = _config[7]; // Use default numi config
		flash_time_end   = _config[8]; // Use default numi config
	}
	std::cout << "flash_time_start:\t" << flash_time_start << "   flash_time_end:\t" << flash_time_end << std::endl;
	
	tolerance = _config[9];
	shwr_nue_tolerance = _config[10];
	trk_nue_tolerance = _config[11];
	shwr_hit_threshold = _config[12];
	shwr_hit_threshold_collection = _config[13];
	tolerance_open_angle_min = _config[14];
	tolerance_open_angle_max = _config[15];
	tolerance_dedx_min = _config[16];
	tolerance_dedx_max = _config[17];
	dist_tolerance = _config[18];
	pfp_hits_length_tolerance = _config[19];
	ratio_tolerance = _config[20];
	detector_variations = _config[21];

	//*************************** SAME PLOT *************************************
	// if bool true just run this function
	if (mode == "same") PlotVar = true; 
	else PlotVar = false; 

	if (PlotVar) {
		PlotVariatons(f_var_out); 
		return; 
	}
	//*************************** POT Scaling *************************************
	std::cout << "=================================================\n" << std::endl;
	
	double CV_POT =  GetPOT("/uboone/data/users/kmistry/work/NueXSection_Outputs/detector_variations/filter_BNBCV.root");
	double POT_Scaling =  CV_POT / GetPOT(_file1);
	std::cout << "POT Scaling:\t" << POT_Scaling << std::endl;
	std::cout << "=================================================\n" << std::endl;
	
	//*********** Open Variation File and get branches ******************************
	// Check if the outfile opened successfully
	if ( f_var_out->IsOpen() ) std::cout << "Variation File opened successfully\n" << std::endl;
	
	// Get Variation file
	TFile* inFile = new TFile(_file1);
	
	// Get TPC Obj Information from file
	TTree* TPCObjTree = (TTree*) inFile->Get("AnalyzeTPCO/tree");	
	std::vector<xsecAna::TPCObjectContainer> * tpc_object_container_v = nullptr;
	TPCObjTree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

	// Num events in the tree
	const int tree_total_entries = TPCObjTree->GetEntries();
	std::cout << "Total Events: " << tree_total_entries << std::endl;

	// Get the largest Flash Vector of Vector	
	std::vector<std::vector<double>> largest_flash_v_v = GetLargestFlashVector(inFile, flash_time_start, flash_time_end);

	//************* Get list of flashes that pass flash in time and flash pe cut ************************
	std::vector<bool> flash_cuts_pass_vec;
	FlashinTime_FlashPE(inFile, flash_time_start, flash_time_end, flash_cuts_pass_vec );
	//***************************************************************************************************

	// MCTruth Counting Tree
	TTree * mctruth_counter_tree = (TTree*)inFile->Get("AnalyzeTPCO/mctruth_counter_tree"); 
	mctruth_counter_tree->SetBranchAddress("fMCNuEnergy", &mc_nu_energy);
	mctruth_counter_tree->SetBranchAddress("fMCNuID", &mc_nu_id);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxX", &mc_nu_vtx_x);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxY", &mc_nu_vtx_y);
	mctruth_counter_tree->SetBranchAddress("fMCNuVtxZ", &mc_nu_vtx_z);
	mctruth_counter_tree->SetBranchAddress("fMCNuDirX", &mc_nu_dir_x);
	mctruth_counter_tree->SetBranchAddress("fMCNuDirY", &mc_nu_dir_y);
	mctruth_counter_tree->SetBranchAddress("fMCNuDirZ", &mc_nu_dir_z);
	mctruth_counter_tree->SetBranchAddress("fMCNumParticles", &mc_nu_num_particles);
	mctruth_counter_tree->SetBranchAddress("fMCNumChargedParticles", &mc_nu_num_charged_particles);
	mctruth_counter_tree->SetBranchAddress("fMCEleDirX", &mc_ele_dir_x);
	mctruth_counter_tree->SetBranchAddress("fMCEleDirY", &mc_ele_dir_y);
	mctruth_counter_tree->SetBranchAddress("fMCEleDirZ", &mc_ele_dir_z);
	mctruth_counter_tree->SetBranchAddress("fMCEleEnergy", &mc_ele_energy);
	mctruth_counter_tree->SetBranchAddress("fMCEleMomentum", &mc_ele_momentum);
	mctruth_counter_tree->SetBranchAddress("has_pi0", &has_pi0);
	mctruth_counter_tree->SetBranchAddress("fMCNuTime", &mc_nu_time);

	// Define the FV
	std::vector<double> fv_boundary_v = {_x1, _x2, _y1, _y2, _z1, _z2};

	std::cout << "Total Events (MC): " << mctruth_counter_tree->GetEntries() << std::endl;

	// ----------------------
	//      Event loop
	// ----------------------
	std::cout << "Starting Eventloop..." << std::endl;
	for(int event = 0; event < tree_total_entries; event++){
		
		TPCObjTree->GetEntry(event);
		mctruth_counter_tree->GetEntry(event);

		int n_tpc_obj = tpc_object_container_v->size();
		if (event < 10) std::cout << "Total TPC Obj: " << n_tpc_obj << std::endl;

		// --------------- MC Counters ---------------
		if(mc_nu_id == 1) mc_nue_cc_counter++;
		if(mc_nu_id == 3) mc_nue_nc_counter++;
		if(mc_nu_id == 5) mc_nue_cc_counter_bar++;
		if(mc_nu_id == 7) mc_nue_nc_counter_bar++;
		
		// Filter for truth Nues
		// if (mc_nu_id == 2 || mc_nu_id == 4 || mc_nu_id == 6 || mc_nu_id == 8) continue;

		// --------------- Flash Information ---------------
		std::vector<double> largest_flash_v 	= largest_flash_v_v.at(event); // Vec with the largest flash
		largest_flash_y 	= largest_flash_v.at(0);
		largest_flash_z 	= largest_flash_v.at(1);
		largest_flash_time 	= largest_flash_v.at(2);
		largest_flash_pe 	= largest_flash_v.at(3);

		// Only fill for matched events (where each entry is non-zero)
		if (largest_flash_y != 0 && largest_flash_z !=0 && largest_flash_time != 0 && largest_flash_pe != 0){
			h_largest_flash_y	->Fill(largest_flash_y);
			h_largest_flash_z	->Fill(largest_flash_z);
			h_largest_flash_time->Fill(largest_flash_time);
			h_largest_flash_pe	->Fill(largest_flash_pe);
		}

		//************************ Apply Flash in time and Flash PE cut *************************************
		if (flash_cuts_pass_vec.at(event) == false) continue;
		//***************************************************************************************************
		
		// if (bool_sig) sig_counter++;
		// -------------------------------------------------

		// Loop over TPCObj
		for (int i = 0; i < n_tpc_obj; i++){
			auto const tpc_obj = tpc_object_container_v->at(i);

			// Check to see in truth the vertex was inside the tpc FV
			const bool true_in_tpc = in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v);

			// Classify the event  
			std::pair<std::string, int> tpc_classification = TPCO_Classifier(tpc_obj, true_in_tpc, has_pi0);

			// Checks the classification for Signal
			bool bool_sig{false}; 
			for (unsigned int k=0; k < signal_modes.size(); k++){
				if (signal_modes.at(k).find(tpc_classification.first) != std::string::npos) {
					bool_sig = true;
				}
			}

			// if (bool_sig) nue_cc_counter++;

			// TPC Obj vars
			tpc_obj_vtx_x = tpc_obj.pfpVtxX();
			tpc_obj_vtx_y = tpc_obj.pfpVtxY();
			tpc_obj_vtx_z = tpc_obj.pfpVtxZ();
		   
			const int tpc_obj_mode = tpc_obj.Mode(); 
			n_pfp = tpc_obj.NumPFParticles();
			const int leading_shower_index = GetLeadingShowerIndex(n_pfp, n_tpc_obj, tpc_obj);
			bool use_secondary_shower_cut{false};
			const int secondary_shower_index = GetSecondaryShowerIndex(n_pfp, n_tpc_obj, tpc_obj);
			if (secondary_shower_index != -999) use_secondary_shower_cut = true;

			// Other histograms
			n_tracks = 0; n_showers = 0; n_pfp_50Hits = 0; n_tracks_50Hits = 0; n_showers_50Hits = 0;
			GetNumber_Track_Shower(n_pfp, n_tpc_obj, tpc_obj, n_showers, n_tracks, n_pfp_50Hits, n_tracks_50Hits, n_showers_50Hits );// Get the number of tracks and showers
			h_n_pfp			  ->Fill(n_tracks + n_showers);
			h_n_pfp_50Hits	  ->Fill(n_pfp_50Hits);
			h_n_tracks		  ->Fill(n_tracks);
			h_n_tracks_50Hits ->Fill(n_tracks_50Hits);
			h_n_showers		  ->Fill(n_showers);
			h_n_showers_50Hits->Fill(n_showers_50Hits);

			const double Flash_TPCObj_Dist = Flash_TPCObj_vtx_Dist(tpc_obj_vtx_y, tpc_obj_vtx_z, largest_flash_y, largest_flash_z);
			h_Flash_TPCObj_Dist->Fill(Flash_TPCObj_Dist);


			//****************************** Apply Pandora Reco Nue cut *****************************************
			bool bool_HasNue = HasNue(tpc_obj, n_pfp );
			if ( bool_HasNue == false ) continue;
			//***************************************************************************************************
			
			//****************************** Apply In FV cut ****************************************************
			bool bool_inFV = in_fv(tpc_obj_vtx_x, tpc_obj_vtx_y, tpc_obj_vtx_z, fv_boundary_v);
			if ( bool_inFV == false ) continue;
			//***************************************************************************************************
			if (bool_sig) nue_cc_counter++;
			

			 

			// Loop over the Par Objects
			for (int j = 0; j < n_pfp ; j++){

				auto const pfp_obj = tpc_obj.GetParticle(j);

				// PFP vars                
				const std::string mc_origin = pfp_obj.Origin();
				const int  pfp_pdg          = pfp_obj.PFParticlePdgCode();
				const int  num_pfp_hits     = pfp_obj.NumPFPHits();
				const int  mc_parent_pdg    = pfp_obj.MCParentPdg();
				const int  n_pfp_hits_w     = pfp_obj.NumPFPHitsW(); // Collection plane hits
				
				const double pfp_length     = pfp_obj.pfpLength();
				const double pfp_open_angle = pfp_obj.pfpOpenAngle();
				const double leading_dedx   = pfp_obj.PfpdEdx().at(2);//just the collection plane!

				const double pfp_vtx_x =  pfp_obj.pfpVtxX();
				const double pfp_vtx_y =  pfp_obj.pfpVtxY();
				const double pfp_vtx_z =  pfp_obj.pfpVtxZ();

				const double pfp_dir_x = pfp_obj.pfpDirX();
				const double pfp_dir_y = pfp_obj.pfpDirY();
				const double pfp_dir_z = pfp_obj.pfpDirZ();

				const double pfp_end_x = (pfp_obj.pfpVtxX() + (pfp_length * pfp_dir_x));
				const double pfp_end_y = (pfp_obj.pfpVtxY() + (pfp_length * pfp_dir_y));
				const double pfp_end_z = (pfp_obj.pfpVtxZ() + (pfp_length * pfp_dir_z));

				std::vector<double> pfp_start_vtx {pfp_vtx_x, pfp_vtx_y, pfp_vtx_z};
				std::vector<double> pfp_end_vtx {pfp_end_x, pfp_end_y, pfp_end_z};

				const double pfp_Nu_vtx_Dist =  pfp_vtx_distance(tpc_obj_vtx_x, tpc_obj_vtx_y, tpc_obj_vtx_z, pfp_vtx_x, pfp_vtx_y, pfp_vtx_z);

				

				// Electron (Shower like)
				if ( pfp_pdg == 11  ) {
					

					// Background events
					if (!bool_sig) {

						// std::cout << tpc_classification.first << std::endl;

						bkg_counter++;
						h_total_hits->Fill(num_pfp_hits);

						const double shower_phi = atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
						h_shower_phi->Fill(shower_phi);

						h_shower_Nu_vtx_Dist->Fill(pfp_Nu_vtx_Dist);

						//  ------------ Leading shower ------------
						if (j == leading_shower_index){
							const double leading_shower_phi 	= atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
							const double leading_shower_theta 	= acos(pfp_obj.pfpDirZ()) * 180 / 3.1415;
							
							const double leading_shower_length 	= pfp_obj.pfpLength();
							const double longest_track_length 	= GetLongestTrackLength(n_pfp, n_tpc_obj, tpc_obj);

							h_ldg_shwr_hits			->Fill(num_pfp_hits);
							h_ldg_shwr_hits_WPlane	->Fill(pfp_obj.NumPFPHitsW()); 		// W Plane
							h_ldg_shwr_Open_Angle	->Fill(pfp_obj.pfpOpenAngle() * (180 / 3.1415) );
							h_ldg_shwr_dEdx_WPlane	->Fill(pfp_obj.PfpdEdx().at(2) ); 	// W Plane
							h_ldg_shwr_HitPerLen	->Fill(num_pfp_hits / pfp_obj.pfpLength() );
							h_ldg_shwr_Phi			->Fill(leading_shower_phi);
							h_ldg_shwr_Theta		->Fill(leading_shower_theta);
							h_ldg_shwr_CTheta		->Fill(cos(leading_shower_theta * 3.1414 / 180.));
							h_long_Track_ldg_shwr	->Fill(longest_track_length / leading_shower_length);

							// Vertex Information - require a shower so fill once  when leading shower
							h_tpc_obj_vtx_x->Fill(tpc_obj_vtx_x);
							h_tpc_obj_vtx_y->Fill(tpc_obj_vtx_y);
							h_tpc_obj_vtx_z->Fill(tpc_obj_vtx_z);
						}
					
					}
					
				}
				// Track like
				if ( pfp_pdg == 13) {
					// Background events
					if (!bool_sig) {
						const double track_phi 	= atan2(pfp_obj.pfpDirY(), pfp_obj.pfpDirX()) * 180 / 3.1415;
						h_track_phi->Fill(track_phi);
						h_track_Nu_vtx_Dist->Fill(pfp_Nu_vtx_Dist);
					}
				}

			} // END LOOP PAR OBJ

		} // END LOOP TPCO

	} // END EVENT LOOP
	std::cout << "Finished Eventloop..." << std::endl;

	
	std::cout << "--------------- MC Truth COUNTERS -----------------" << std::endl;
	std::cout << "MC Nue CC Counter      --- " << mc_nue_cc_counter << std::endl;
	std::cout << "MC Nue NC Counter      --- " << mc_nue_nc_counter << std::endl;
	std::cout << "MC Nue CC Counter Bar  --- " << mc_nue_cc_counter_bar << std::endl;
	std::cout << "MC Nue NC Counter Bar  --- " << mc_nue_nc_counter_bar << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;
	std::cout << "--------------- Reco COUNTERS ---------------------" << std::endl;
	std::cout << "(Requiring shower like and (nue/nuebar or numu/numubar))" << std::endl;
	std::cout << "RECO Nue CC Counter (Tot) --- " << nue_cc_counter << std::endl;
	std::cout << "RECO Signal Counter       --- " << sig_counter << std::endl;
	std::cout << "RECO Background Counter   --- " << bkg_counter << std::endl;
	std::cout << "---------------------------------------------------" << std::endl;
	
	// ----------------------
	//    Save to a file
	// ----------------------
	TDirectory* savedir = gDirectory; //  Create the directory
	TDirectory* subdir;

	// Get the directory for the file	
	std::string dirname = _file1;
	
	// Get the variation name by stripping the input file name -- HARDCODED
	// Required format "files/filter_<NAME>.root"
	std::string file_path = "/uboone/data/users/kmistry/work/NueXSection_Outputs/detector_variations/";
	int file_path_size = file_path.length();
	dirname.erase(0,file_path_size+7); // - "files/filter_"
	dirname.erase(dirname.end()-5, dirname.end()); // - ".root"

	std::cout << "dirname:\t" << dirname << std::endl;

	f_var_out->cd();

	// If directory does not exist then make it
	savedir = (TDirectory*)f_var_out->Get(dirname.c_str());
	if (savedir == NULL ) {
		savedir = gDirectory;
		std::cout << dirname << " directory does not exist, creating..." << std::endl;
		subdir = savedir->mkdir(dirname.c_str() ) ;
	}
	else {
		std::cout << dirname <<" directory exists, overwriting..." << std::endl;
		subdir = savedir;
	}
	
	subdir->cd();

	// ----------------------
	//    Draw Histograms
	// ----------------------
	DrawTH1D(h_total_hits, POT_Scaling);
	DrawTH1D(h_ldg_shwr_hits, POT_Scaling);
	DrawTH1D(h_ldg_shwr_hits_WPlane, POT_Scaling);
	DrawTH1D(h_ldg_shwr_Open_Angle, POT_Scaling);
	DrawTH1D(h_ldg_shwr_dEdx_WPlane, POT_Scaling);
	DrawTH1D(h_ldg_shwr_HitPerLen, POT_Scaling);
	DrawTH1D(h_ldg_shwr_Phi, POT_Scaling);
	DrawTH1D(h_ldg_shwr_Theta, POT_Scaling);
	DrawTH1D(h_ldg_shwr_CTheta, POT_Scaling);
	DrawTH1D(h_long_Track_ldg_shwr, POT_Scaling);
	DrawTH1D(h_tpc_obj_vtx_x, POT_Scaling);
	DrawTH1D(h_tpc_obj_vtx_y, POT_Scaling);	
	DrawTH1D(h_tpc_obj_vtx_z, POT_Scaling);	
	DrawTH1D(h_n_pfp, POT_Scaling);
	DrawTH1D(h_n_pfp_50Hits, POT_Scaling);
	DrawTH1D(h_n_tracks, POT_Scaling);
	DrawTH1D(h_n_tracks_50Hits, POT_Scaling);
	DrawTH1D(h_n_showers, POT_Scaling);
	DrawTH1D(h_n_showers_50Hits, POT_Scaling);
	DrawTH1D(h_track_phi, POT_Scaling);
	DrawTH1D(h_shower_phi, POT_Scaling);
	DrawTH1D(h_largest_flash_y, POT_Scaling);
	DrawTH1D(h_largest_flash_z, POT_Scaling);
	DrawTH1D(h_largest_flash_time, POT_Scaling);
	DrawTH1D(h_largest_flash_pe, POT_Scaling);
	DrawTH1D(h_Flash_TPCObj_Dist, POT_Scaling);
	DrawTH1D(h_shower_Nu_vtx_Dist, POT_Scaling);
	DrawTH1D(h_track_Nu_vtx_Dist, POT_Scaling);

	// ----------------------
	//   Write and close
	//   the TFile to new/updated
	//   directory
	// ----------------------
	h_total_hits->Write("",TObject::kOverwrite);
	h_ldg_shwr_hits->Write("",TObject::kOverwrite);
	h_ldg_shwr_hits_WPlane-> Write("",TObject::kOverwrite);
	h_ldg_shwr_Open_Angle->Write("", TObject::kOverwrite);
	h_ldg_shwr_dEdx_WPlane->Write("", TObject::kOverwrite);
	h_ldg_shwr_HitPerLen->Write("", TObject::kOverwrite);
	h_ldg_shwr_Phi->Write("", TObject::kOverwrite);
	h_ldg_shwr_Theta->Write("", TObject::kOverwrite);
	h_ldg_shwr_CTheta->Write("", TObject::kOverwrite);
	h_long_Track_ldg_shwr->Write("", TObject::kOverwrite);
	h_tpc_obj_vtx_x->Write("", TObject::kOverwrite);
	h_tpc_obj_vtx_y->Write("", TObject::kOverwrite);
	h_tpc_obj_vtx_z->Write("", TObject::kOverwrite);
	h_n_pfp->Write("", TObject::kOverwrite);
	h_n_pfp_50Hits->Write("", TObject::kOverwrite);
	h_n_tracks->Write("", TObject::kOverwrite);
	h_n_tracks_50Hits->Write("", TObject::kOverwrite);
	h_n_showers->Write("", TObject::kOverwrite);
	h_n_showers_50Hits->Write("", TObject::kOverwrite);
	h_track_phi->Write("", TObject::kOverwrite);
	h_shower_phi->Write("", TObject::kOverwrite);
	h_largest_flash_y->Write("", TObject::kOverwrite);
	h_largest_flash_z->Write("", TObject::kOverwrite);
	h_largest_flash_time->Write("", TObject::kOverwrite);
	h_largest_flash_pe->Write("", TObject::kOverwrite);
	h_Flash_TPCObj_Dist->Write("", TObject::kOverwrite);
	h_shower_Nu_vtx_Dist->Write("", TObject::kOverwrite);
	h_track_Nu_vtx_Dist->Write("", TObject::kOverwrite);
	
	f_var_out->Close(); 

} // END MAIN

