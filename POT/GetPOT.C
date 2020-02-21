#include <iostream>
#include <fstream>

//Root Includes
#include "TFile.h"
#include "TTree.h"

/*
Script to get the POT from a pandora ntuple file

Usage: root -l -q -b 'GetPOT.C("/path/to/file","mc/data")'

if MC is given, then it will look for the pot branch in the file, if data is given,
then Zarko's POT counting script will be called after generating a run subrun file list

*/


int GetPOT(const char *_file1, std::string type){

    bool debug = false;

    std::cout << "File: " << _file1 << std::endl;

    if (type == "mc"){
        // First we need to open the root file
        TFile * f = new TFile(_file1);
        if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }

        // Get the tree
        TTree * mytree = (TTree*)f->Get("nuselection/SubRun");

        if (mytree == NULL){
            std::cout << "help can't get the branch so exiting..." << std::endl;
            gSystem->Exit(1);
        }

        float pot_sum = 0;
        float pot;
        mytree->SetBranchAddress("pot", &pot);

        // Loop over tree and get POT
        for (int i = 0; i < mytree->GetEntries(); i++) {
            mytree->GetEntry(i);
            if (debug) std::cout << pot << std::endl;
            pot_sum = pot_sum + pot;
        }

        std::cout << "Total POT: " << pot_sum  << std::endl;
        return 0;

    }
    // Data file so do zarko's POT counting script
    else {

        int run, subrun;

        std::ofstream run_subrun_file;
        run_subrun_file.open("run_subrun_list_data.txt");

        // First we need to open the root file
        TFile * f = new TFile(_file1);
        if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }

        // Get the tree
        TTree * mytree = (TTree*)f->Get("nuselection/SubRun");

        if (mytree == NULL){
            std::cout << "help can't get the branch so exiting..." << std::endl;
            gSystem->Exit(1);
        }

        mytree->SetBranchAddress("run", &run);
        mytree->SetBranchAddress("subRun", &subrun);

        for (int i = 0; i < mytree->GetEntries(); i++){
            mytree->GetEntry(i);
            run_subrun_file << run << " " << subrun << '\n';
        }
        
        run_subrun_file.close();

        gSystem->Exec("/uboone/app/users/zarko/getDataInfo.py -v2 --format-numi --prescale --run-subrun-list run_subrun_list_data.txt"); 
        return 0;

    }
    
}