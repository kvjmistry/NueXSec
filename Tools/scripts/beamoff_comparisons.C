#include "functions.h"

/* 

Script to check the beamoff normalisation for a couple of variables

*/


void beamoff_comparisons(){

    TFile *f_ext, *f_beamoff;

    TH1D *h_ext, *h_beamoff;

    GetFile(f_ext, "/uboone/app/users/kmistry/MCC9_uboonecode_v08_00_00_33/srcs/ubana/ubana/NueXSec/Analysis/files/nuexsec_ext_run1_full.root");
    // GetFile(f_beamoff, "/uboone/app/users/kmistry/MCC9_uboonecode_v08_00_00_33/srcs/ubana/ubana/NueXSec/Analysis/files/nuexsec_ext_run1_filter_absolute.root");
    GetFile(f_beamoff, "/uboone/app/users/kmistry/MCC9_uboonecode_v08_00_00_33/srcs/ubana/ubana/NueXSec/Analysis/files/nuexsec_ext_run1_filter_should.root");


    GetHist(f_ext, h_ext, "Flash/h_flash_time_EXT");
    GetHist(f_beamoff, h_beamoff, "Flash/h_flash_time_EXT");

    h_ext->Rebin(5);
    h_beamoff->Rebin(5);

    // h_beamoff->Scale(2987202.860000 / 42562.0); // abolute sample
    h_beamoff->Scale(2987202.860000 / 55913.0);  // should sample

    gStyle->SetOptStat(0);

    TLegend *leg_stack = new TLegend();
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);

    leg_stack->AddEntry(h_ext, "EXT",   "f");
    // leg_stack->AddEntry(h_beamoff, "Beam On No Beam",   "l"); // absolute samole
    leg_stack->AddEntry(h_beamoff, "Beam On Should Sample",   "l"); // should sample

    TCanvas *c = new TCanvas();

    double y_maximum = std::max(h_ext->GetMaximum(), h_beamoff->GetMaximum());

    h_ext->GetYaxis()->SetRangeUser(0, y_maximum+3000);
    // h_beamoff->GetYaxis()->SetRangeUser(0, 25000);

    h_beamoff->SetLineColor(kBlack);

    h_ext      ->SetFillColor(41);
    h_ext      ->SetFillStyle(3345);
    h_ext->Draw("hist,E");
    h_beamoff->Draw("same,E");
    leg_stack->Draw();

    const char* test = "empty";

    std::string test2 = test;

    // std::string hi = "files/" + test;


}