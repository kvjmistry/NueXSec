#include "../include/selection_cuts.h"
// -----------------------------------------------------------------------------
void selection_cuts::SetFlashVariables(std::vector<double> largest_flash_v){
        largest_flash_y 	= largest_flash_v.at(0);
        largest_flash_z 	= largest_flash_v.at(1);
        largest_flash_time 	= largest_flash_v.at(2);
        largest_flash_pe 	= largest_flash_v.at(3);
}
// -----------------------------------------------------------------------------
void SetTPCObjVariables(xsecAna::TPCObjectContainer tpc_obj){



}
// -----------------------------------------------------------------------------