// ROOT script to plot the flash info 

void plot_flashinfo(){

    gStyle->SetOptStat(0);

    // First get the files
    TFile * f_Overlay = TFile::Open("/uboone/data/users/kmistry/work/MCC9/FlashValidation/flash_validation_run1_overlay.root");
    TFile * f_Dirt    = TFile::Open("/uboone/data/users/kmistry/work/MCC9/FlashValidation/flash_validation_run1_dirt_overlay.root");
    TFile * f_Data    = TFile::Open("/uboone/data/users/kmistry/work/MCC9/FlashValidation/flash_validation_run1_data.root");
    TFile * f_EXT     = TFile::Open("/uboone/data/users/kmistry/work/MCC9/FlashValidation/flash_validation_run1_ext.root");

    // Now get the histograms
    TH1D *h_flash_overlay = (TH1D*) f_Overlay->Get("FlashValidate/h_flash_time");
    TH1D *h_flash_dirt    = (TH1D*) f_Dirt   ->Get("FlashValidate/h_flash_time");
    TH1D *h_flash_data    = (TH1D*) f_Data   ->Get("FlashValidate/h_flash_time");
    TH1D *h_flash_ext     = (TH1D*) f_EXT    ->Get("FlashValidate/h_flash_time");

    TH1D *h_flash_pe_overlay = (TH1D*) f_Overlay->Get("FlashValidate/h_flash_pe");
    TH1D *h_flash_pe_dirt    = (TH1D*) f_Dirt   ->Get("FlashValidate/h_flash_pe");
    TH1D *h_flash_pe_data    = (TH1D*) f_Data   ->Get("FlashValidate/h_flash_pe");
    TH1D *h_flash_pe_ext     = (TH1D*) f_EXT    ->Get("FlashValidate/h_flash_pe");

    if (h_flash_overlay == NULL) std::cout << "Error couldnt get the overlay histogram" << std::endl;
    if (h_flash_dirt == NULL)    std::cout << "Error couldnt get the dirt histogram" << std::endl;
    if (h_flash_data == NULL)    std::cout << "Error couldnt get the data histogram" << std::endl;
    if (h_flash_ext == NULL)     std::cout << "Error couldnt get the ext histogram" << std::endl;

    // Define POT scalings and HW trigger amounts for each sample
    double MC_POT    = 3.56691e+20;
    double Dirt_POT  = 2.85818e+20;
    double Data_POT  = 3.918e+19;
    double Data_trig = 1009015.0;
    double EXT_trig  = 851387.820000;

    // double scale_factor = 0.8;

    // Scale the histograms
    h_flash_overlay ->Scale(scale_factor * Data_POT / MC_POT);
    h_flash_dirt    ->Scale(scale_factor * Data_POT / Dirt_POT);
    h_flash_ext     ->Scale(scale_factor * Data_trig / EXT_trig);

    h_flash_pe_overlay ->Scale(Data_POT / MC_POT);
    h_flash_pe_dirt    ->Scale(Data_POT / Dirt_POT);
    h_flash_pe_ext     ->Scale(Data_trig / EXT_trig);

    // Make the stack
    THStack * h_stack    = new THStack();
    THStack * h_pe_stack = new THStack();

    h_flash_overlay->SetFillColor(30);
    h_flash_dirt   ->SetFillColor(38);
    h_flash_ext    ->SetFillColor(41);
    h_flash_ext    ->SetFillStyle(3345);
        
    h_flash_pe_overlay->SetFillColor(30);
    h_flash_pe_dirt   ->SetFillColor(38);
    h_flash_pe_ext    ->SetFillColor(41);
    h_flash_pe_ext    ->SetFillStyle(3345);

    h_stack->Add(h_flash_ext);
    h_stack->Add(h_flash_overlay);
    h_stack->Add(h_flash_dirt);

    h_pe_stack->Add(h_flash_pe_ext);
    h_pe_stack->Add(h_flash_pe_overlay);
    h_pe_stack->Add(h_flash_pe_dirt);

    //h_stack->GetYaxis()->SetTitleFont(45);
    //h_stack->GetYaxis()->SetTitleSize(18);
    //h_stack->GetYaxis()->SetTitleOffset(1.30);
    //h_stack->GetXaxis()->SetLabelOffset(10);
    
    TLegend *leg_stack = new TLegend(0.70, 0.60, 0.90, 0.90);
    //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header	
    leg_stack->AddEntry(h_flash_data,     "Data",    "l");
    leg_stack->AddEntry(h_flash_ext,      "EXT",     "f");
    leg_stack->AddEntry(h_flash_overlay,  "Overlay", "f");
    leg_stack->AddEntry(h_flash_dirt,     "Dirt",    "f");

    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);
    
    TCanvas* c       = new TCanvas();
    h_flash_data->GetYaxis()->SetRangeUser(0,2100);
    h_flash_data->Draw("same PE");
    h_stack->Draw("hist,same");
    h_flash_data->Draw("same PE");
    leg_stack ->Draw();
    c->Print("Flash_Time.pdf");

    TCanvas* c2       = new TCanvas();
    h_pe_stack->SetTitle("Flash PE; Flash PE [PE]; Entries");
    h_pe_stack->Draw("hist");
    h_flash_pe_data->Draw("same PE");
    leg_stack ->Draw();
    c2->Print("Flash_PE.pdf");

}
