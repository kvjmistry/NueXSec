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
double _Run1_MC_POT    = 3.37458e+20;      // Run1 MC POT v33

// double _Run1_MC_POT    = 1.26178e+20;       // Run1 MC POT v33 --elena's sample
double _Run1_Dirt_POT  = 6.43292e+20;       // Run 1 dirt POT


// double _Run1_Data_POT  = 4.822e+19;         // Run 1 Data POT (tortgt_wcut)
// double _Run1_Data_trig = 1278305.0;         // Run 1 data HW Triggers (EA9CNT_wcut) 
double _Run1_Data_POT  =  3.932e+19 ;          // Run 1 Data POT (tortgt) beam good and bad file
double _Run1_Data_trig = 1052383.0;            // Run 1 data HW Triggers (EA9CNT)  beam good and bad file

double _Run1_EXT_trig  = 2987202.860000;    // Run 1 Number of EXT HW Triggers ( EXT_NUMIwin_FEMBeamTriggerAlgo )

// Fiducial Volume
const double _x1 = 10;
const double _x2 = 246.35;
const double _y1 = -106.5;
const double _y2 = 106.5;
const double _z1 = 20;
const double _z2 = 986.8;

#endif
