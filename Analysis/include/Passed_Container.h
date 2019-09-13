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
        Passed_Container()=default;
        
        // Fiducial volume
        bool in_fv{true};

        // Flash PE
        bool flash_pe{true};

        // In time flash
        bool flash_intime{true};

        // Vertex to flash
        bool vtx_to_flash{true};

        // Distance between pfp shower and nue object
        bool shwr_nue_dist{true};
        
        // Distance between pfp track and nue object
        bool trk_nue_dist{true};

        // Hit threshold for at least one shower
        bool shwr_hit_threshold{true};

        // Hit threshold for at least one shower on collection plane
        bool shwr_hit_threshold_collection{true};

        // Tolerance for leading shower open angle
        bool shwr_open_angle{true};

        // Tolerance for dedx of leading shower
        bool shwr_dedx{true};

        // Tolerance for distance from the reco nue vtx for TPCO w/ >3 showers
        bool dist_nue_vtx{true};

        // Tolerance for hits/length
        bool pfp_hits_length{true};

        // Tolerance for longest track length / leading shower length
        bool longest_trk_leading_shwr_length{true};
    
};

#endif