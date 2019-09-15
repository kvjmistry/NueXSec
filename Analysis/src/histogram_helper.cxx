#include "../include/histogram_helper.h"

// -----------------------------------------------------------------------------
histogram_helper::~histogram_helper() { 
    
    // Make sure the file is closed
    f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void histogram_helper::MakeDirectory(std::string type){
        
    f_nuexsec->cd();

    TDirectory *top_dir; // e.g MC, Data, EXT
    bool bool_dir;       // Check if directory exists already
    TString type_tstr = type;
   
    // Create the top directory
    bool_dir = _utility_instance.GetDirectory(f_nuexsec, top_dir, type_tstr );
    if (!bool_dir) top_dir = f_nuexsec->mkdir(type.c_str());
    
    // Make the the top dir the current directory
    top_dir->cd();

    // Create subdirectory for cut type
    TDirectory *dir_plot_types[plot_types.size()];
    
    // Create a new subdirectory for each cut
    const Int_t ncuts = cut_dirs.size();
    TDirectory *dir_cut[ncuts];

    // Create a new subdirectory for each classification
    TDirectory *dir_classification[classification_dirs.size()];
    
    // Loop over the plot types ------------------------------------------------
    for (unsigned int k = 0; k < plot_types.size(); k++) {
        
        // Get the directory 
        bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_plot_types[k] ,Form("%s/%s", type.c_str(), plot_types.at(k).c_str()) );

        // Make the directory
        if (!bool_dir) dir_plot_types[k] = top_dir->mkdir(plot_types.at(k).c_str());

        dir_plot_types[k]->cd();

        // If we have stacked histograms, we make plots by cut 
        if (plot_types.at(k) == "Stack"){
            
            // Loop over the cuts ----------------------------------------------
            for (int i = 0; i < ncuts; i++) {
               
                // Get the directory 
                bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_cut[i] ,Form("%s/%s/%s", type.c_str(), plot_types.at(k).c_str(), cut_dirs.at(i).c_str()));

                // Make the directory
                if (!bool_dir) dir_cut[i] = dir_plot_types[k]->mkdir(cut_dirs.at(i).c_str());
                dir_cut[i]->cd();
                
                // Loop over the classifications -------------------------------
                for (unsigned int j = 0; j < classification_dirs.size(); j++){
                
                    // Get the directory 
                    bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_classification[j] ,Form("%s/%s/%s/%s", type.c_str(), plot_types.at(k).c_str(), cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()));

                    // Make the directory
                    if (!bool_dir) dir_classification[j] = dir_cut[i]->mkdir(classification_dirs.at(j).c_str());
                    dir_classification[j]->cd();
                    top_dir->cd();    // change current directory to top

                } // End loop over the classifications -------------------------
                
                top_dir->cd();    // change current directory to top
                
            } // End loop over the cuts ----------------------------------------
       
        }
         
        top_dir->cd();    // change current directory to top
    
    } // End loop over plot types ----------------------------------------------

    top_dir->Write("",TObject::kOverwrite);
   
}
// -----------------------------------------------------------------------------
void histogram_helper::Initialise(){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject("nuexsec.root") ) {
        f_nuexsec = new TFile("nuexsec.root", "UPDATE");
    }

    MakeDirectory("MC");
    MakeDirectory("Data");
    MakeDirectory("EXT");
    MakeDirectory("Dirt");
}
// -----------------------------------------------------------------------------