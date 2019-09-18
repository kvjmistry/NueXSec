#ifndef HISTOGRAM_PLOTTER_h
#define HISTOGRAM_PLOTTER_h

#include "histogram_helper.h"

// Class for getting histograms which have been written to file and plot them 
// in a pretty way
// Inherits all tools from the histogram helper
class histogram_plotter : public histogram_helper{

    public:
    // Default constructor
    histogram_plotter(){};
    
    // Destructor 
    ~histogram_plotter(); 


    // -------------------------------------------------------------------------
    // Initialise histograms
    void InitHistograms();
    // -------------------------------------------------------------------------
    // Template to making a stacked histogram. Sets the relavent variables up
    // void SetStack(std::vector<TH1D*> hist, bool area_norm,  bool logy, const char* x_axis_name,
    //                                  double data_scale_factor, double y_scale_factor, double intime_scale_factor, double dirt_scale_factor, 
    //                                  const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, const char* print_name );
    // -------------------------------------------------------------------------
    // Returns the chi squared
    // std::vector <double> Chi2Calc(TH1D * h_mc_ext, TH1D * h_data, const bool area_norm, const double return_norm);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    private:

    std::vector<THStack*> h_reco_vtx_x_stack;
    std::vector<THStack*> h_reco_vtx_y_stack;
    std::vector<THStack*> h_reco_vtx_z_stack;

   

}; // End Class Histogram Plotter

#endif
