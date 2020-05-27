#include <iostream>
#include <fstream>

//Root Includes
#include "TFile.h"
#include "TTree.h"

/*
Script to make the plot of run number and POT normalisation

Usage: root -l -q -b 'RunPeriod_Norm.C("<merged histogram file from selection>)'
*/

// Cut directory names -- needs to be kept up to date with the selection
std::vector<std::string> cut_dirs = {
        "Unselected",     // Unselected
        "SoftwareTrig",   // Software Trigger
        "Op_Filter_PE",   // Common Optical Filter PE
        "Op_Filter_Veto",// Common Optical Filter Michel Veto
        "Slice_ID",       // Slice ID
        "e_candidate",    // Electron Candidate
        "In_FV",          // In FV
        "Topo_Score",     // Topological Score
        "Cosmic_IP",      // Cosmic Inpact Parameter
        "Cluster_Frac",   // Cluster Fraction 
        "Contained_Frac", // Slice Contained Fraction
        "Shower_Score",   // Track Score
        "Michel_Rej",     // Michel Rejection
        "ShrHits",        // Shower Hits
        "HitRatio",       // Ratio of shr hits and slice hits
        "Moliere_Avg",    // Shower Moliere Average
        "ShrVtxDistance", // Shower to vertex distance
        "dEdx_y",         // dEdx y plane
        };
// -----------------------------------------------------------------------------
void RunPeriod_Norm(const char *_file1){

    std::string run_period = "run1";

    gStyle->SetOptStat(0);

    std::cout << "File: " << _file1 << std::endl;

    
    // First we need to open the root file
    TFile * f = new TFile(_file1);
    if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; }

    // Get the histograms
    std::vector<TH1D*> h_run_data_EA9CNT;
    std::vector<TH1D*> h_run_data_pot_v;
    std::vector<TH1D*> h_run_ext;

    h_run_data_EA9CNT.resize(cut_dirs.size());
    h_run_data_pot_v.resize(cut_dirs.size());
    h_run_ext.resize(cut_dirs.size());

    // Get the histograms of run number as a function of number of events
    for (unsigned int l=0; l < h_run_data_EA9CNT.size(); l++){
        h_run_data_EA9CNT.at(l) = (TH1D*) f->Get( Form("Stack/%s/data/h_reco_run_number_%s_data",cut_dirs.at(l).c_str() ,cut_dirs.at(l).c_str()));
        h_run_data_pot_v.at(l)  = (TH1D*) h_run_data_EA9CNT.at(l)->Clone(Form("h_run_data_clone_%s", cut_dirs.at(l).c_str())  );
        h_run_ext.at(l)         = (TH1D*) f->Get( Form("Stack/%s/ext/h_reco_run_number_%s_ext", cut_dirs.at(l).c_str() ,cut_dirs.at(l).c_str()));
    }
    
    // Now we need to make the run number as a function of POT, EA9CNT and EXT Triggers
    int run, subrun;
    double pot, EA9CNT, ext_trig;
    std::string line;

    TH1D* h_run_ext_trig  = new TH1D("h_run_ext_trig",  "; Run Number; EXT Triggers",                     270, 4500, 18000);
    TH1D* h_run_data_trig = new TH1D("h_run_data_EA9CNT_trig", "; Run Number; EA9CNT (Data HW Triggers)", 270, 4500, 18000);
    TH1D* h_run_data_pot  = new TH1D("h_run_data_EA9CNT_pot",  "; Run Number; tortgt_wcut (Data POT)",    270, 4500, 18000);

    // EXT Triggers
    std::ifstream myfile ("filelists/beamoff_trig.txt");
    if (myfile.is_open()) {
        
        // Loop over lines in file
        while ( getline (myfile,line) ) {

            std::istringstream ss(line);
            ss >> run >> subrun >> pot >> EA9CNT >> ext_trig;

            h_run_ext_trig->Fill(run, ext_trig);
            
        }
    }
    else std::cout << "Unable to open file, bad things are going to happen..." << std::endl; 

    myfile.close();

    // Data POT and triggers
    std::ifstream myfile2 ("filelists/beamon_pot_trig.txt");
    if (myfile2.is_open()) {
        
        // Loop over lines in file
        while ( getline (myfile2,line) ) {

            std::istringstream ss(line);
            ss >> run >> subrun >> pot >> EA9CNT >> ext_trig;

            if (EA9CNT != 0 || pot != 0 ){
                h_run_data_trig->Fill(run, EA9CNT);
                h_run_data_pot->Fill(run, pot);
            }
        }
    }
    else std::cout << "Unable to open file, bad things are going to happen..." << std::endl;

    myfile2.close();


    // Now got everything, take the ratios!
    for (unsigned int l=0; l < h_run_data_EA9CNT.size(); l++){
        h_run_data_EA9CNT.at(l)->Divide(h_run_data_trig);
        h_run_data_pot_v.at(l)->Divide(h_run_data_pot);
        h_run_ext.at(l)->Divide(h_run_ext_trig);
    }

    // Make the EXT plot
    TCanvas *c = new TCanvas();

    if (run_period == "run1") h_run_ext.at(0)->GetXaxis()->SetRangeUser(4500, 8000);
    else h_run_ext.at(0)->GetXaxis()->SetRangeUser(16880, 18000);

    h_run_ext.at(0)->SetLineWidth(2);
    h_run_ext.at(0)->GetYaxis()->SetTitle("Events / EXT Triggers");
    h_run_ext.at(0)->SetLineColor(kBlack);
    h_run_ext.at(0)->Draw("E");

    h_run_ext.at(3)->SetLineWidth(2);
    h_run_ext.at(3)->SetLineColor(kRed+2);
    h_run_ext.at(3)->Draw("E,same");

    h_run_ext.at(4)->SetLineWidth(2);
    h_run_ext.at(4)->SetLineColor(kBlue+2);
    h_run_ext.at(4)->Draw("E,same");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_run_ext.at(0), "No Cuts", "l");
    leg->AddEntry(h_run_ext.at(3), "Optical Filter", "l");
    leg->AddEntry(h_run_ext.at(4), "Slice ID", "l");
    
    leg->Draw();

    c->Print("plots/ext_events_per_trig.pdf");

    // Make the EA9CNT plot
    TCanvas *c2 = new TCanvas();

    if (run_period == "run1") h_run_data_EA9CNT.at(0)->GetXaxis()->SetRangeUser(4500, 8000);
    else h_run_data_EA9CNT.at(0)->GetXaxis()->SetRangeUser(16880, 18000);

    h_run_data_EA9CNT.at(0)->SetLineWidth(2);
    h_run_data_EA9CNT.at(0)->GetYaxis()->SetTitle("Events / EA9CNT (Data Triggers)");
    h_run_data_EA9CNT.at(0)->SetLineColor(kBlack);
    h_run_data_EA9CNT.at(0)->Draw("E");

    h_run_data_EA9CNT.at(3)->SetLineWidth(2);
    h_run_data_EA9CNT.at(3)->SetLineColor(kRed+2);
    h_run_data_EA9CNT.at(3)->Draw("E,same");

    h_run_data_EA9CNT.at(4)->SetLineWidth(2);
    h_run_data_EA9CNT.at(4)->SetLineColor(kBlue+2);
    h_run_data_EA9CNT.at(4)->Draw("E,same");

    leg->Draw();

    c2->Print("plots/data_events_per_trig.pdf");

    TCanvas *c3 = new TCanvas();

    if (run_period == "run1") h_run_data_pot_v.at(0)->GetXaxis()->SetRangeUser(4500, 8000);
    else h_run_data_pot_v.at(0)->GetXaxis()->SetRangeUser(16880, 18000);

    h_run_data_pot_v.at(0)->SetLineWidth(2);
    h_run_data_pot_v.at(0)->GetYaxis()->SetTitle("Events / tortgt_wcut (Data POT)");
    h_run_data_pot_v.at(0)->SetLineColor(kBlack);
    h_run_data_pot_v.at(0)->Draw("E");

    h_run_data_pot_v.at(3)->SetLineWidth(2);
    h_run_data_pot_v.at(3)->SetLineColor(kRed+2);
    h_run_data_pot_v.at(3)->Draw("E,same");

    h_run_data_pot_v.at(4)->SetLineWidth(2);
    h_run_data_pot_v.at(4)->SetLineColor(kBlue+2);
    h_run_data_pot_v.at(4)->Draw("E,same");

    leg->Draw();

    c3->Print("plots/data_events_per_pot.pdf");

    TCanvas *c4 = new TCanvas();
    h_run_data_pot->Draw();





 
}
// -----------------------------------------------------------------------------