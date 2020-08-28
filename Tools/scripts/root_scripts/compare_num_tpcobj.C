// Script to plot the stacked version of the tpc objects

// enums for cut dirs
enum enum_cut_dirs {
            k_unselected,        // Unselected 
            k_in_time,
            k_flash_pe,
            k_reco_nue,
            k_in_fv,
            k_vtx_to_flash,
            k_shr_to_vtx,
            k_track_to_vtx,
            k_hit_thresh,
            k_hit_thresh_y,
            k_open_angle,
            k_dedx,
            k_shr_vtx_dist_gt_1shr,
            k_hit_per_lenth,
            k_trk_shr_lengh,
            k_trk_containment,
            k_cuts_MAX
            }; 

enum enum_type {
            k_mc,
            k_dirt,
            k_data,
            k_ext
            }; 




// Cut directory names
std::vector<std::string> cut_dirs = {
        "Unselected",     // Unselected
        "In_Time",
        "Flash_PE",
        "Reco_Nue",
        "In_FV",
        "Vtx_to_Flash",
        "Shr_to_Vtx",
        "Track_to_Vtx",
        "Hit_Thresh",
        "Hit_Thresh_Y",
        "Open_Angle",
        "dEdx",
        "Shr_Vtx_Dist_gt_1shr",
        "Hit_per_Lenth",
        "Trk_Shr_Lengh",
        "Trk_Containment"
        };




void compare_num_tpcobj(){

    // Load in the Tfiles
    TFile *f_mc   = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_mc_mcc8.root");
    TFile *f_dirt = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_dirt_mcc8.root");
    TFile *f_data = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_data_mcc8.root");
    TFile *f_ext  = TFile::Open("/Users/kvjmistry/Documents/work/MCC8/NueXSec/files/f_ext_mcc8.root");

    double mc_scale_factor     =  0.1301;
    double dirt_scale_factor   =  0.16411; //0.45
    double intime_scale_factor =  1.0154; // 0.98


    // Define the histograms
    std::vector<std::vector<TH1D*>> h_obj;

    h_obj.resize(4); // mc/data/dirt/ext

    for (unsigned int j = 0; j < h_obj.size(); j++){
        h_obj.at(j).resize(k_cuts_MAX);
    }

    // Now lets get the histgrams
    f_mc->cd();
    for (int k = 0; k < k_cuts_MAX; k++){
        h_obj.at(k_mc).at(k) = (TH1D*)f_mc->Get(Form("h_%s_ntpc", cut_dirs.at(k).c_str()));
        h_obj.at(k_mc).at(k)->Scale(mc_scale_factor);
        h_obj.at(k_mc).at(k)->SetFillColor(30);

    }

    f_data->cd();
    for (int k = 0; k < k_cuts_MAX; k++){
        h_obj.at(k_data).at(k) =(TH1D*) f_data->Get(Form("h_%s_ntpc", cut_dirs.at(k).c_str()));
        h_obj.at(k_data).at(k)->SetMarkerStyle(20);
        h_obj.at(k_data).at(k)->SetMarkerSize(0.5);
    }

    f_dirt->cd();
    for (int k = 0; k < k_cuts_MAX; k++){
        h_obj.at(k_dirt).at(k) = (TH1D*)f_dirt->Get(Form("h_%s_ntpc", cut_dirs.at(k).c_str()));
        h_obj.at(k_dirt).at(k)->Scale(dirt_scale_factor);
        h_obj.at(k_dirt).at(k)->SetFillColor(2);
        h_obj.at(k_dirt).at(k)->SetFillStyle(3354);
    }

    f_ext->cd();
    for (int k = 0; k < k_cuts_MAX; k++){
        h_obj.at(k_ext).at(k) =(TH1D*) f_ext->Get(Form("h_%s_ntpc", cut_dirs.at(k).c_str()));
        h_obj.at(k_ext).at(k)->Scale(intime_scale_factor);
        h_obj.at(k_ext).at(k)->SetFillColor(41);
        h_obj.at(k_ext).at(k)->SetFillStyle(3345);
    }

    TLegend *leg_stack = new TLegend(0.5, 0.6, 0.8, 0.8);
    leg_stack->SetBorderSize(0);
    leg_stack->SetFillStyle(0);
    leg_stack->AddEntry(h_obj.at(k_data).at(0), "Data", "lep");
    leg_stack->AddEntry(h_obj.at(k_dirt).at(0), "Dirt", "f");
    leg_stack->AddEntry(h_obj.at(k_ext).at(0), "EXT", "f");
    leg_stack->AddEntry(h_obj.at(k_mc).at(0), "MC", "f");

    
    for (int k = 0; k < k_cuts_MAX; k++){
        TCanvas *c = new TCanvas("c", "c", 500, 500);

        THStack *h_stack = new THStack();

        h_stack->Add(h_obj.at(k_mc).at(k));
        h_stack->Add(h_obj.at(k_ext).at(k));
        h_stack->Add(h_obj.at(k_dirt).at(k));
        h_stack->SetMinimum(0);

        h_obj.at(k_data).at(k)->SetTitle(";Num TPC Objects; Entries");
        
        h_obj.at(k_data).at(k)->Draw("E,same");
        h_stack->Draw("hist,same");
        h_obj.at(k_data).at(k)->Draw("E,same");
        leg_stack->Draw();


        c->Print(Form("plots/%i_h_%s_ntpc.pdf", k, cut_dirs.at(k).c_str()));
        delete c;
    }




}