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

    // Scale factors (veverything is scaled to data)
    double mc_scale_factor, intime_scale_factor, dirt_scale_factor;
    
    bool weight_tune = true; // Apply genie tune weight
    bool weight_ppfx = true; // Apply ppfx cv weight
    bool area_norm   = false;  // Decide if to area normalse the histograms


    // -------------------------------------------------------------------------
    // Main function call to control this class
    void MakeHistograms(const char * hist_file_name, const char *run_period, int weight_cfg, bool _area_norm, utility _utility, const char* variation);
    // -------------------------------------------------------------------------
    // Initalise the file input
    void Initalise(const char *hist_file_name, const char* _run_period, double _mc_scale_factor, double _intime_scale_factor, double _dirt_scale_factor, int weight_cfg);
    // -------------------------------------------------------------------------
    // Function to make a stacked histogram and save as a pdf
    void MakeStack(std::string hist_name, std::string cut_name, bool area_norm, bool logy, double y_scale_factor, const char* x_axis_name,
                                     const double leg_x1, const double leg_x2, const double leg_y1, double Data_POT, const double leg_y2,
                                     const char* print_name, bool override_data_mc_comparison, std::string plotmode, bool plotvar,
                                     const char * variation, const char *run_period, bool centerxaxis );
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
    // Draw area norm text 
    void Draw_Area_Norm(TCanvas* c);
    // -------------------------------------------------------------------------
    // Add the labels for weights
    void Draw_WeightLabels(TCanvas* c);
    // -------------------------------------------------------------------------
    // Draw Variation Mode
    void Draw_VarMode(TCanvas* c, const char* variation);
    // -------------------------------------------------------------------------
    // Set the TPad Options
    void SetTPadOptions(TPad * topPad, TPad * bottomPad );
    // -------------------------------------------------------------------------
    // Increase the label size of the 1D histogram
    void IncreaseLabelSize(TH1D* h);
    // -------------------------------------------------------------------------
    // Increase the label size of the 2D histogram
    void IncreaseLabelSize(TH2D* h);
    // -------------------------------------------------------------------------
    // Function to get all the histograms from the file
    bool GetHistograms(std::vector<TH1D*> &hist, std::string hist_name, std::string cut_name, std::string plotmode, bool &found_data, bool &found_ext, bool &found_dirt);
    // -------------------------------------------------------------------------
    // Function to set the colours of the stacked histogram
    void SetFillColours(std::vector<TH1D*> &hist, std::string plotmode, bool found_data, bool found_dirt, bool found_ext, unsigned int k_plot_data, unsigned int k_plot_ext, unsigned int k_plot_dirt);
    // -------------------------------------------------------------------------
    // Function to set the legend of the stacked histograms
    void SetLegend(std::vector<TH1D*> hist, TLegend *leg_stack, std::vector<double> hist_integrals, bool found_data, bool found_dirt, bool found_ext, unsigned int k_plot_data, unsigned int k_plot_ext, unsigned int k_plot_dirt, std::string plotmode );
    // -------------------------------------------------------------------------
    // Call make stacked histograms
    void CallMakeStack(const char *run_period, int cut_index, double Data_POT, const char *variation);
    // -------------------------------------------------------------------------
    void MakeFlashPlot(double Data_POT, const char* print_name, std::string histname);
    // -------------------------------------------------------------------------
    // Flash plot for on beam minus off beam
    void MakeFlashPlotOMO(double Data_POT, const char* print_name, std::string histname);
    // -------------------------------------------------------------------------
    // Make efficiency plots
    void MakeEfficiencyPlot(const char* print_name, const char *run_period);
    // -------------------------------------------------------------------------
    // The efficiency made by each cut
    void MakeEfficiencyPlotByCut(const char *run_period);
    // -------------------------------------------------------------------------
    // True Neutino energy for nues broken down by genie interaction type
    void MakeInteractionPlot(const char* print_name);
    // -------------------------------------------------------------------------
    // Plot the 2D signal vs Background Plots
    void Plot2D_Signal_Background(const char* print_name, const char* histname);
    // -------------------------------------------------------------------------
    // Create another directory in the plots folder
    void CreateDirectory(std::string folder, const char *run_period);
    // -------------------------------------------------------------------------
    // Script to get the 1D histograms and save them as a PDF
    void Save1DHists(const char* print_name, const char* histname);
    // -------------------------------------------------------------------------
    // Script to get the 2D histograms and save them as a PDF
    void Save2DHists(const char* print_name, const char* histname);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------


    private:


}; // End Class Histogram Plotter

#endif
