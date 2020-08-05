#ifndef PASSEDCONTAINER_H
#define PASSEDCONTAINER_H

#include <vector>
#include <string>

#include "Utility.h"

/* 
Class to hold information to whether a specific cut has passed the selection
Main purpose is to be a container.

We basically use this container to write a filelist of passed events so we can look at them
in the event display.

Default is to set true and change to false when cut fails

*/

// Passed Container Class
class PassedContainer {
    public:

        Utility _util;

        PassedContainer();
        
        std::vector<bool> cut_v;
        
};

#endif