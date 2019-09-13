#ifndef PASSED_CONTAINER_H
#define PASSED_CONTAINER_H

#include <vector>
#include <string>

/* 
Class to hold information to whether a specific cut has passed the selection
Main purpose is to be a container.

Default is to set true and change to false when cut fails

*/

// Passed Container Class
class Passed_Container {
    public:
        Passed_Container(){
            cut_v.resize(k_cuts_MAX, true);
        };
        
        enum cuts{
            k_in_fv,                           // Fiducial volume
            k_flash_pe,                        // Flash PE
            k_flash_intime,                    // In time flash
            k_vtx_to_flash,                    // Vertex to flash
            k_shwr_nue_dist,                   // Distance between pfp shower and nue object
            k_trk_nue_dist,                    // Distance between pfp track and nue object
            k_shwr_hit_threshold,              // Hit threshold for at least one shower
            k_shwr_hit_threshold_collection,   // Hit threshold for at least one shower on collection plane
            k_shwr_open_angle,                 // Tolerance for leading shower open angle
            k_shwr_dedx,                       // Tolerance for dedx of leading shower
            k_dist_nue_vtx,                    // Tolerance for distance from the reco nue vtx for TPCO w/ >3 showers
            k_pfp_hits_length,                 // Tolerance for hits/length
            k_longest_trk_leading_shwr_length, // Tolerance for longest track length / leading shower length
            k_cuts_MAX
        };

        std::vector<bool> cut_v;

};

#endif