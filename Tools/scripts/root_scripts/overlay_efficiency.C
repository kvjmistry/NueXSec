
// Root script to plot the efficiency and purity as a function of cuts for each run period to compare them

#include "functions.h"

// Cut directory names -- harcoded, so need to be updated
std::vector<std::string> cut_dirs = {
        "Unselected",     // Unselected
        "SoftwareTrig",   // Software Trigger
        "Slice_ID",       // Slice ID
        "e_candidate",    // Electron Candidate
        "In_FV",          // In FV
        "Topo_Score",     // Topological Score
        "Contained_Frac", // Slice Contained Fraction
        "Shower_Score",   // Track Score
        "Michel_Rej",     // Michel Rejection
        "ShrHits",        // Shower Hits
        "HitRatio",       // Ratio of shr hits and slice hits
        "Moliere_Avg",    // Shower Moliere Average
        "ShrVtxDist_dEdx_y",
        "dEdx_y_no_tracks",  // dEdx y plane no tracks
        // "ShrVtxDistance", // Shower to vertex distance
        // "dEdx_y",         // dEdx y plane
        };

// enums for cut dirs
enum enum_cut_dirs {
            k_unselected,        // Unselected 
            k_swtrig,            // Software Trigger
            k_slice_id,          // Slice ID
            k_e_candidate,       // Electron Candidate
            k_in_fv,             // Reco Nu Vtx (SC Corr) In the FV 
            k_topo_score,        // Topo Score
            k_contained_frac,    // Slice Contained Fraction
            k_shower_score,      // Shower Score
            k_michel_rej,        // Michel Rejection
            k_shr_hits,          // Shower Hits
            k_hit_ratio,         // Ratio of shr hits and slice hits
            k_shr_moliere_avg,   // Shower Moliere Average
            k_vtx_dist_dedx,     //  2D cut for shower to vertex distance and dEdx. Only applied for > 1 track
            k_dEdx_y_no_tracks,  // dEdx y plane when there is no tracks
            // k_shr_distance,      // Shower to reco nu vertex distance
            // k_dEdx_y,            // dEdx y plane
            k_cuts_MAX
            }; 

void populate_efficiency_vec(TH1D* h_eff, TH1D *h_pur, const char* input_file){

    TTree *mc_tree;

    TFile *f;

    // The uglyest hardcoded monstosity known to the universe. SORT this out KRISH...
    f = TFile::Open(input_file);

    std::vector<double> efficiency_v; // efficiency vector
    std::vector<double> purity_v    ; // purity vector

    double efficiency, purity;

    GetTree(f, mc_tree, "mc_eff_tree");
    mc_tree->SetBranchAddress("efficiency", &efficiency);
    mc_tree->SetBranchAddress("purity",     &purity);

    int num_entries = mc_tree->GetEntries();

    // Fill the efficiency and purity vectors
    for (int y=0; y < num_entries; y++){
        mc_tree->GetEntry(y); 
        efficiency_v.push_back(efficiency);
        purity_v.push_back(purity);
    }

    for (unsigned int k=0; k < efficiency_v.size();k++){
        h_eff ->Fill(cut_dirs.at(k).c_str(), efficiency_v.at(k));
        h_pur ->Fill(cut_dirs.at(k).c_str(), purity_v.at(k));
        h_eff->SetBinError(k+1, 0);
        h_pur->SetBinError(k+1, 0);
    }

    h_eff->GetYaxis()->SetRangeUser(0,1.1);
    h_eff->SetStats(kFALSE);
    h_eff->SetMarkerStyle(20);
    h_eff->SetMarkerSize(0.5);
    h_eff->SetLineWidth(2);

    h_pur->SetLineColor(kRed+2);
    h_pur->SetStats(kFALSE);
    h_pur->SetMarkerStyle(20);
    h_pur->SetMarkerSize(0.5);
    h_pur->SetLineWidth(2);


}


void overlay_efficiency(){
    
    TH1D* h_eff_run1 = new TH1D("h_efficiency_run1", "", k_cuts_MAX, 0, k_cuts_MAX);
    TH1D* h_pur_run1 = new TH1D("h_purity_run1",     "", k_cuts_MAX, 0, k_cuts_MAX);

    TH1D* h_eff_run3 = new TH1D("h_efficiency_run3", "", k_cuts_MAX, 0, k_cuts_MAX);
    TH1D* h_pur_run3 = new TH1D("h_purity_run3",     "", k_cuts_MAX, 0, k_cuts_MAX);

    populate_efficiency_vec(h_eff_run1, h_pur_run1, "../../../Analysis/files/trees/nuexsec_selected_tree_mc_run1.root");
    populate_efficiency_vec(h_eff_run3, h_pur_run3, "../../../Analysis/files/trees/nuexsec_selected_tree_mc_run3.root");
    
    
    TCanvas *c = new TCanvas();
    c->SetGridy();

    TLegend *leg_stack = new TLegend(0.7, 0.9, 0.90, 0.7);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);
    leg_stack->AddEntry(h_eff_run1, "Efficiency Run1","lp");
    leg_stack->AddEntry(h_eff_run3, "Efficiency Run3","lp");
    leg_stack->AddEntry(h_pur_run1, "Purity Run1",    "lp");
    leg_stack->AddEntry(h_pur_run3, "Purity Run3",    "lp");

    h_eff_run1->Draw("LP");
    h_pur_run1->Draw("LP,same");

    h_eff_run3->SetLineColor(kAzure+5);
    h_eff_run3->Draw("LP,same");

    h_pur_run3->SetLineColor(kViolet-5);
    h_pur_run3->Draw("LP,same");
    
    leg_stack->Draw();

    // Draw vertical lines to help the eye
    TLine *line;
    for (unsigned int l=1; l < k_cuts_MAX+1; l++){
        line  = new TLine( h_eff_run1->GetBinCenter(l) ,   0 , h_eff_run1->GetBinCenter(l)  ,  1.1);
        line->SetLineColor(12);
        line->SetLineStyle(kDotted);
        line->Draw();
    }

    h_eff_run1->GetXaxis()->SetTickLength(0.00);
    h_pur_run1->GetXaxis()->SetTickLength(0.00);

    c->Print("../../../Analysis/plots/efficiency_run1_run3_comparison.pdf");



}