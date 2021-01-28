#include <iostream>
#include <fstream>
#include <string.h>
#include "TString.h"

//Root Includes
#include "TFile.h"
#include "TTree.h"

/*
Script to get the POT from a pandora ntuple file

Usage: root -l -q -b 'GetPOT.C("/path/to/file","mc/data/ext")'

if MC is given, then it will look for the pot branch in the file, if data is given,
then Zarko's POT counting script will be called after generating a run subrun file list

*/


void GetPOT(const char *_file1, std::string type, int run_number){

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

	// ----- check if the POT number in config.txt matched with pot_sum

	std::cout << std::endl;
	std::cout << "Checking the POT value: " << std::endl;

	ifstream config_file("../Analysis/config.txt"); // file with the POT values
	bool same_pot = false;

	std::string line; // saves the each line from config.txt
	
	while(getline(config_file,line)){
		
		if(line[0]=='R'){ // skips comment lines

			// splitting the line and taking the first arg as a TString and the second one as a float
			TString line_split(line);
			TObjArray *substrings = line_split.Tokenize(" ");
			float val = ((TObjString*)substrings->At(1))->GetString().Atof();
			TString label = ((TObjString*)substrings->At(0))->GetString();

			// calculate the ratio between new and old value, if =1 they are the same
			double ratio_pot = pot_sum/val;
			string ratio_pot_str = std::to_string(ratio_pot);

			// checks the initial label and if the ratio is one
			if(label == "Run1_MC_POT_CV" && ratio_pot_str == "1.000000") { // going back to string gives me a precision of 6 decimals
				std::cout << "The POT value of your file match the one in config.txt, continue running the code!" << std::endl;
				std::cout << std::endl;
			}

			// if POT values don't match, stop the code
			else if(label == "Run1_MC_POT_CV" && ratio_pot_str != "1.000000") {
				std::cout << "########## Warning message !!! ##########" << std::endl;
				std::cout << "POT value from your input file: " << pot_sum << std::endl;
				std::cout << "POT value saved in config.txt:  " << val << std::endl;
				std::cout << "The POT value of your file does not match the one in config.txt." << std::endl;
				std::cout << "#########################################" << std::endl;
				std::cout << std::endl;
				throw std::invalid_argument("The POT values don't match, please update config.txt");
			}
		}
	}

    }
    // Data file so do zarko's POT counting script
    else {

        int run, subrun;

        std::cout << "     "  << std::endl;
        system("if test -f \"run_subrun_list_data.txt\"; then echo \"the run_subrun file list exists, removing as a safeguard...\"; fi");
        system("if test -f \"run_subrun_list_data.txt\"; then rm run_subrun_list_data.txt; fi");
        std::cout << "     "  << std::endl;

        std::ofstream run_subrun_file;
        run_subrun_file.open("run_subrun_list_data.txt");

        // First we need to open the root file
        TFile * f = new TFile(_file1);
        if(!f->IsOpen()) {std::cout << "Could not open file!" << std::endl; return 1; }

        // Get the tree
        TTree * mytree = (TTree*)f->Get("nuselection/SubRun");
        if (mytree == NULL) mytree = (TTree*)f->Get("FlashValidate/pottree");

        if (mytree == NULL){
            std::cout << "help can't get the branch so exiting..." << std::endl;
            gSystem->Exit(1);
        }

        mytree->SetBranchAddress("run", &run);
        mytree->SetBranchAddress("subRun", &subrun);
        //mytree->SetBranchAddress("subrun", &subrun);

        for (int i = 0; i < mytree->GetEntries(); i++){
            mytree->GetEntry(i);

            //if (run < 16880) continue;
            //if (run > 16880) continue;
            //if (run >6748.22) continue;
            //
            //if (run < 7013) continue;
            if (run == 16228) continue;

            run_subrun_file << run << " " << subrun << '\n';
        }
        
        run_subrun_file.close();

        if (type == "ext") {
            std::cout << "/uboone/app/users/zarko/getDataInfo.py -v3 --format-numi --prescale --run-subrun-list run_subrun_list_data.txt" << std::endl;
            gSystem->Exec("/uboone/app/users/zarko/getDataInfo.py -v3 --format-numi --prescale --run-subrun-list run_subrun_list_data.txt"); 
        }
        // Use Pawels updated version for on beam, use  --slip for slipstacking info
        else {
            std::cout << "/uboone/app/users/guzowski/slip_stacking/getDataInfo.py -v3 --format-numi --prescale --run-subrun-list run_subrun_list_data.txt --slip" << std::endl;
            gSystem->Exec("/uboone/app/users/guzowski/slip_stacking/getDataInfo.py -v3 --format-numi --prescale --run-subrun-list run_subrun_list_data.txt --slip"); 

        } 
    
    
        gSystem->Exec("rm run_subrun_list_data.txt");

    }
    
}
