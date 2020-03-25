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

    std::string run_period;

    double mc_scale_factor, intime_scale_factor, dirt_scale_factor;

    // -------------------------------------------------------------------------
    // Initalise the file input
    void Initalise(const char *hist_file_name, const char* _run_period, double _mc_scale_factor, double _intime_scale_factor, double _dirt_scale_factor);
    // -------------------------------------------------------------------------
    // Function to make a stacked histogram and save as a pdf
    void MakeStack(std::string hist_name, std::string cut_name, bool area_norm, bool logy, double y_scale_factor, const char* x_axis_name,
                                     const double leg_x1, const double leg_x2, const double leg_y1, double Data_POT, const double leg_y2, const char* print_name );
    // -------------------------------------------------------------------------
    // Calculates the chi2
    // Returns reduced chi2, num mc+ext scaled to data POT, num data, num degrees of freedom, p value in vector respectively
    std::vector<double> Chi2Calc(TH1D * h_mc_ext, TH1D * h_data, const bool area_norm, const double return_norm);
    // -------------------------------------------------------------------------
    // Draw the run period on the plot
    void Draw_Run_Period(TCanvas* c);
    // -------------------------------------------------------------------------
    // Draw the Ratio to MC ratio
    void Draw_Data_MC_Ratio(TCanvas* c, double ratio);
    // -------------------------------------------------------------------------
    // Draw the data POT
    void Draw_Data_POT(TCanvas* c, double pot);
    // -------------------------------------------------------------------------
    // Call make stacked histograms
    void CallMakeStack(const char *run_period, int cut_index, double Data_POT);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    private:


}; // End Class Histogram Plotter

#endif
