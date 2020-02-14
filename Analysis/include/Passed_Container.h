#ifndef PASSED_CONTAINER_H
#define PASSED_CONTAINER_H

#include <vector>
#include <string>

/* 
Class to hold information to whether a specific cut has passed the selection
Main purpose is to be a container.

Will also hold counters

Default is to set true and change to false when cut fails

*/

// Passed Container Class
class Passed_Container {
    public:
        Passed_Container(){
            cut_v.resize(k_cuts_MAX, false);
            
        };
        
        enum cuts{
            k_in_fv,                           // Fiducial volume
            k_cuts_MAX
        };

        std::vector<bool> cut_v;
        
};

#endif