#ifndef PASSED_CONTAINER_H
#define PASSED_CONTAINER_H

#include <vector>
#include <string>

#include "utility.h"

/* 
Class to hold information to whether a specific cut has passed the selection
Main purpose is to be a container.

Will also hold counters

Default is to set true and change to false when cut fails

*/

// Passed Container Class
class Passed_Container {
    public:

        utility _util;

        Passed_Container();
        
        std::vector<bool> cut_v;
        
};

#endif