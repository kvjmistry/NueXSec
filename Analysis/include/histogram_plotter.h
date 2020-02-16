#ifndef HISTOGRAM_PLOTTER_h
#define HISTOGRAM_PLOTTER_h

#include "utility.h"

// Class for getting histograms which have been written to file and plot them 
// in a pretty way
// Inherits all tools from the histogram helper
class histogram_plotter{

    public:
    
    
    // Destructor 
    ~histogram_plotter(); 

    TFile *f_nuexsec;

    utility _util;

    // -------------------------------------------------------------------------
    // Default constructor
    void Initalise(int type);
    // -------------------------------------------------------------------------
    // Initialise histograms
    void InitHistograms();
    // -------------------------------------------------------------------------
    // Function to make a stacked histogram and save as a pdf
    void MakeStack(std::string hist_name, std::string cut_name, bool area_norm,  bool logy, const char* x_axis_name,
                                     double data_scale_factor, double y_scale_factor, double intime_scale_factor, double dirt_scale_factor, 
                                     const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, const char* print_name );
    // -------------------------------------------------------------------------
    // Calculates the chi2
    // Returns reduced chi2, num mc+ext scaled to data POT, num data, num degrees of freedom, p value in vector respectively
    std::vector<double> Chi2Calc(TH1D * h_mc_ext, TH1D * h_data, const bool area_norm, const double return_norm);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    private:


}; // End Class Histogram Plotter

#endif
