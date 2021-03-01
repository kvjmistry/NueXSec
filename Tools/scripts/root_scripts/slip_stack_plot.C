


void SetHistProperties(TCanvas *c, TH1D* h){

    h->GetXaxis()->SetLabelSize(0.03);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.11);

    // Draw the histogram
    h->SetLineWidth(0);
    h->SetStats(kFALSE);

}

void SetHistProperties(TCanvas *c, TH2D* h){

    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(1.8);
    h->GetZaxis()->SetLabelSize(0.05);
    h->GetZaxis()->SetTitleSize(0.05);
    c->SetLeftMargin(0.2);
    c->SetRightMargin(0.2);
    c->SetBottomMargin(0.13);
    h->SetMarkerSize(1.8);

}


void slip_stack_plot(){

    gStyle->SetOptStat(0);

    TFile *f = TFile::Open("test.root", "RECREATE");

    TTree* rrstree = new TTree("runsubtree","Run SubRun Tree");
    TTree* tree4p6 = new TTree("4p6","4 + 6 TTree");
    TTree* tree6p6 = new TTree("6p6","6 + 6 TTree");
    TTree* treeother = new TTree("other","Other TTree");
    rrstree->ReadFile("files/run_subrun_list_data.txt","run/I:subrun/I");
    tree4p6->ReadFile("files/beamon_slip_4p6.txt","pot4p6/D");
    tree6p6->ReadFile("files/beamon_slip_6p6.txt","pot6p6/D");
    treeother->ReadFile("files/beamon_slip_other.txt","other/D");
    
    int _run, _subrun;
    double _pot4p6, _pot6p6, _other;

    int run, subrun;
    double pot4p6, pot6p6, other;

    rrstree->SetBranchAddress("run", &_run);
    rrstree->SetBranchAddress("subrun", &_subrun);
    tree4p6->SetBranchAddress("pot4p6", &_pot4p6);
    tree6p6->SetBranchAddress("pot6p6", &_pot6p6);
    treeother->SetBranchAddress("other", &_other);
    


    TTree* slipstacktree = new TTree("slipstacktree","TTree with Slip Stack info");
    slipstacktree->Branch("run", &run);
    slipstacktree->Branch("subrun", &subrun);
    slipstacktree->Branch("pot4p6", &pot4p6);
    slipstacktree->Branch("pot6p6", &pot6p6);
    slipstacktree->Branch("other", &other);

    
    int tree_total_entries = rrstree->GetEntries();

    // Event loop
    for (int ievent = 0; ievent < tree_total_entries; ievent++){

        // Alert the user
        if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;

        // Get the entry in the tree
        rrstree->GetEntry(ievent); 
        tree4p6->GetEntry(ievent); 
        tree6p6->GetEntry(ievent); 
        treeother->GetEntry(ievent); 

        run = _run;
        subrun = _subrun;
        pot4p6 = _pot4p6;
        pot6p6 = _pot6p6;
        other = _other;

        slipstacktree->Fill();
    }


    // Now Lets query the TTree and make some plots!
    TCanvas * c = new TCanvas("c", "c", 500, 500);
    
    TH1D* htemp4p6 = new TH1D("htemp4p6",";Run Number; POT ", 50, 4900, 6900);
    slipstacktree->Draw(" run >> htemp4p6", "pot4p6");
    SetHistProperties(c, htemp4p6);
    htemp4p6->SetFillColor(46);
    

    TH1D* htemp6p6 = new TH1D("htemp6p6",";Run Number; POT ", 50, 4900, 6900);
    slipstacktree->Draw(" run >> htemp6p6", "pot6p6");
    SetHistProperties(c, htemp6p6);
    htemp6p6->SetFillColor(30);
    
    THStack *h_stack = new THStack();
    h_stack->SetTitle(";Run Number; POT ");
    
    h_stack->Add(htemp4p6);
    h_stack->Add(htemp6p6);

    h_stack->Draw("hist");
    h_stack->GetXaxis()->SetLabelSize(0.03);
    h_stack->GetXaxis()->SetTitleSize(0.05);
    h_stack->GetYaxis()->SetLabelSize(0.05);
    h_stack->GetYaxis()->SetTitleSize(0.05);

    // Create the legend
    TLegend *leg = new TLegend(0.2, 0.70, 0.50, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(htemp4p6, "4+6 Slip Stacking", "f");
    leg->AddEntry(htemp6p6, "6+6 Slip Stacking", "f");
    leg->Draw();

    TPaveText *pt;
    pt = new TPaveText(0.86, 0.92, 0.86, 0.92, "NDC");
    pt->AddText("Run1");
    pt->SetTextColor(kRed + 2);
    pt->SetTextSize(0.04);
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.04);
    pt->Draw();

    c->Print("plots/run1_slip_stacking.pdf");

    // Now Lets query the TTree and make some plots!
    TCanvas * c2 = new TCanvas("c2", "c2", 500, 500);
    TH2D* htemprun4p6 = new TH2D("htemprun4p6",";Run Number; 4+6 Intensity", 15, 4900, 6900, 15, 0, 8e15);
    slipstacktree->Draw(" pot4p6:run  >> htemprun4p6", "pot4p6", "colz");
    SetHistProperties(c2, htemprun4p6);
    c2->SetTopMargin(0.11);
    htemprun4p6->GetXaxis()->SetLabelSize(0.03);
    gStyle->SetPalette(kBlueGreenYellow);
    c2->Print("plots/run1_4p6_slip_stacking_Intensity.pdf");

    // Now Lets query the TTree and make some plots!
    TCanvas * c3 = new TCanvas("c3", "c3", 500, 500);
    TH2D* htemprun6p6 = new TH2D("htemprun6p6",";Run Number; 6+6 Intensity", 15, 4900, 6900, 15, 0, 8e15);
    slipstacktree->Draw(" pot6p6:run  >> htemprun6p6", "pot6p6", "colz");
    SetHistProperties(c3, htemprun6p6);
    c3->SetTopMargin(0.11);
    htemprun6p6->GetXaxis()->SetLabelSize(0.03);
    gStyle->SetPalette(kBlueGreenYellow);
    c3->Print("plots/run1_6p6_slip_stacking_Intensity.pdf");

    slipstacktree->Write();

    f->Close();
}