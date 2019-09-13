#include "../include/utility.h"

namespace utilityNS {

    // -------------------------------------------------------------------------
    std::vector<double> utility::configure_cuts(
            double _x1,
            double _x2,
            double _y1,
            double _y2,
            double _z1,
            double _z2,
            double flash_pe_threshold,
            double flash_time_start,
            double flash_time_end,
            double tolerance,
            double shwr_nue_tolerance,
            double trk_nue_tolerance,
            double shwr_hit_threshold,
            double shwr_hit_threshold_collection,
            double tolerance_open_angle_min,
            double tolerance_open_angle_max,
            double tolerance_dedx_min,
            double tolerance_dedx_max,
            double dist_tolerance,
            double pfp_hits_length_tolerance,
            double ratio_tolerance,
            bool   do_variations
            ) {
        std::vector<double> config;
        config.resize(22,0);

        config[0]  = _x1;
        config[1]  = _x2;
        config[2]  = _y1;
        config[3]  = _y2;
        config[4]  = _z1;
        config[5]  = _z2;
        config[6]  = flash_pe_threshold;
        config[7]  = flash_time_start;
        config[8]  = flash_time_end;
        config[9]  = tolerance;
        config[10] = shwr_nue_tolerance;
        config[11] = trk_nue_tolerance;
        config[12] = shwr_hit_threshold;
        config[13] = shwr_hit_threshold_collection;
        config[14] = tolerance_open_angle_min;
        config[15] = tolerance_open_angle_max;
        config[16] = tolerance_dedx_min;
        config[17] = tolerance_dedx_max;
        config[18] = dist_tolerance;
        config[19] = pfp_hits_length_tolerance;
        config[20] = ratio_tolerance;
        config[21] = do_variations;

        return config;

    } // End config function
    // -------------------------------------------------------------------------
    bool utility::GetFile(TFile* &f, TString string){
        f = TFile::Open(string);
        
        if (f == NULL) {
            std::cout << "failed to get:\t" << string << "\tThis mode wont be used in the selection" << std::endl;
            return false;
        }
        else {
            return true;
        }
    }
    // -------------------------------------------------------------------------
    void utility::GetTree(TFile* f, TTree* &T, TString string){
        T = (TTree*) f->Get(string);
        if (T == NULL) {
            std::cout << "\nfailed to get:\t" << string << "\tThis tree might not exist in the file, exiting...\n" << std::endl;
            exit(1);
        }
        else {
            return;
        }
    }
    // -------------------------------------------------------------------------
    bool utility::in_fv(double x, double y, double z, std::vector<double> fv_boundary_v) {
    
        // The MicroBooNE Active Volume
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

        // Not in AV
        if(x <= det_x1 + x1 || x >= det_x2 - x2) return false;
        if(y <= det_y1 + y1 || y >= det_y2 - y2) return false;
        if(z <= det_z1 + z1 || z >= det_z2 - z2) return false;
        
        // Passed
        return true;
    }
    // -------------------------------------------------------------------------
    void utility::GetLargestFlashVector(TTree* optical_tree, double flash_time_start, double flash_time_end, int flash_pe_threshold, 
        std::vector<std::vector<double>> &largest_flash_v_v, std::vector<std::vector<int>> &optical_list_pe_v, std::vector<std::vector<double>> &optical_list_flash_time_v){

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
    std::cout << "Total Optical Entries:     " << optical_entries << std::endl;

    int current_event = 0;
    int current_run = 0;
    int last_event = 0;
    int last_run = 0;

    // Contains the entry number for a given OpFlash per event
    std::vector<int>					optical_list_pe;
    
    std::vector<double>					optical_list_flash_center_y; 
    std::vector<std::vector<double> >	optical_list_flash_center_y_v;
    
    std::vector<double>					optical_list_flash_center_z; 
    std::vector<std::vector<double> >	optical_list_flash_center_z_v;
    
    std::vector<double>					optical_list_flash_time;
    
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
    for(unsigned int i = 0; i < optical_list_pe_v.size(); i++) {
        
        bool in_time 				= false;
        bool sufficient_flash 		= false;
        
        double largest_flash = 0.;
        double largest_center_y = 0;
        double largest_center_z = 0;
        double largest_flash_time = 0;
        
        // Cut Variables defined in main.h
        int flash_pe_threshold = 50;

        // Loop through all flashes in event and find largest
        for(unsigned int j = 0; j < optical_list_pe_v.at(i).size(); j++) {
            
            auto const opt_time         = optical_list_flash_time_v.at(i).at(j) ; // shift due to MC and offset
            auto const opt_pe           = optical_list_pe_v.at(i).at(j);
            const double opt_center_y   = optical_list_flash_center_y_v.at(i).at(j);
            const double opt_center_z   = optical_list_flash_center_z_v.at(i).at(j);
            const double opt_flash_time = optical_list_flash_time_v.at(i).at(j); // shift due to MC and offset
            
            // See if flash was in time
            in_time = (opt_time >= flash_time_start && opt_time <= flash_time_end) ? true : false;
            
            // See if flash meets the threshold requirements
            sufficient_flash = (opt_pe >= flash_pe_threshold) ? true : false;
            
            // Flash is both in time and over PE threshold
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

    return;
    }
    // -------------------------------------------------------------------------
} // End namespace utlity
