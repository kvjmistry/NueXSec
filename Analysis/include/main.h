#ifndef MAIN_h
#define MAIN_h

#include "utility.h"
#include "selection.h"

/*

This is the main header

In this header we configure any changable values in the selection. I aim to make 
sure that there is nothing hardcoded in the code further down.

*/

// POT and Trigger 
double _Run1_MC_POT    = 1.66215e+21;      // Run1 MC POT v33

// double _Run1_MC_POT    = 1.26178e+20;       // Run1 MC POT v33 --elena's sample
double _Run1_Dirt_POT  = 3.82938e+20;       // Run 1 dirt POT

// double _Run1_Data_POT  = 4.822e+19;         // Run 1 Data POT (tortgt_wcut)
// double _Run1_Data_trig = 1278305.0;         // Run 1 data HW Triggers (EA9CNT_wcut) 
double _Run1_Data_POT  =  3.932e+19 ;          // Run 1 Data POT (tortgt) beam good and bad file
double _Run1_Data_trig = 1052383.0;            // Run 1 data HW Triggers (EA9CNT)  beam good and bad file

double _Run1_EXT_trig  = 3378479.250000;    // Run 1 Number of EXT HW Triggers ( EXT_NUMIwin_FEMBeamTriggerAlgo )

// POT and Trigger 
double _Run3_MC_POT    = 1.14486e+21;      // Run 3 MC POT v33
double _Run3_Dirt_POT  = 4.32295e+20;      // Run 3 dirt POT
double _Run3_Data_POT  = 4.049e+19;         // Run 3 Data POT (tortgt) beam good and bad file
double _Run3_Data_trig = 827436.0 ;         // Run 3 data HW Triggers (EA9CNT)  beam good and bad file
double _Run3_EXT_trig  = 1558577.475000;   // Run 1 Number of EXT HW Triggers ( EXT_NUMIwin_FEMBeamTriggerAlgo )



// Fiducial Volume
const double _x1 = 10;
const double _x2 = 246.35;
const double _y1 = -106.5;
const double _y2 = 106.5;
const double _z1 = 20;
const double _z2 = 986.8;

#endif
