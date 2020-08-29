#ifndef HISTOGRAMPLOTTER_H
#define HISTOGRAMPLOTTER_H

#include "Utility.h"

// Class for getting histograms which have been written to file and plot them 
// in a pretty way
// Inherits all tools from the histogram helper
class HistogramPlotter{

    public:
    
    
    // Destructor 
    ~HistogramPlotter(); 

    TFile *f_nuexsec;

    Utility _util;

    bool area_norm   = false;  // Decide if to area normalse the histograms


    // -------------------------------------------------------------------------
    // Main function call to control this class
    void MakeHistograms(Utility _utility);
    // -------------------------------------------------------------------------
    // Initalise the file input
    void Initalise();
    // -------------------------------------------------------------------------
    // Function to make a stacked histogram and save as a pdf
    void MakeStack(std::string hist_name, std::string cut_name, bool area_norm, bool logy, double y_scale_factor, const char* x_axis_name,
                                     const double leg_x1, const double leg_x2, const double leg_y1, double Data_POT, const double leg_y2,
                                     const char* print_name, bool override_data_mc_comparison, std::string plotmode, bool plotvar,
                                    bool centerxaxis, bool scale );
    // -------------------------------------------------------------------------
    // Calculates the chi2
    // Returns reduced chi2, num mc+ext scaled to data POT, num data, num degrees of freedom, p value in vector respectively
    std::vector<double> Chi2Calc(TH1D * h_mc_ext, TH1D * h_data, const bool area_norm, const double return_norm);
    // -------------------------------------------------------------------------
    // Draw area norm text 
    void Draw_Area_Norm(TCanvas* c);
    // -------------------------------------------------------------------------
    // Add the labels for weights
    void Draw_WeightLabels(TCanvas* c);
    // -------------------------------------------------------------------------
    // Draw Variation Mode
    void Draw_VarMode(TCanvas* c);
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
    void CallMakeStack(int cut_index, double Data_POT);
    // -------------------------------------------------------------------------
    void MakeFlashPlot(double Data_POT, const char* print_name, std::string histname);
    // -------------------------------------------------------------------------
    // Flash plot for on beam minus off beam
    void MakeFlashPlotOMO(double Data_POT, const char* print_name, std::string histname);
    // -------------------------------------------------------------------------
    // Make efficiency plots
    void MakeEfficiencyPlot(const char* print_name);
    // -------------------------------------------------------------------------
    // The efficiency made by each cut
    void MakeEfficiencyPlotByCut(std::string var, bool mask_title);
    // -------------------------------------------------------------------------
    // True Neutino energy for nues broken down by genie interaction type
    void MakeInteractionPlot(const char* print_name, std::string cut_type, bool scale);
    // -------------------------------------------------------------------------
    // Plot the 2D signal vs Background Plots
    void Plot2D_Signal_Background(const char* print_name, const char* histname);
    // -------------------------------------------------------------------------
    // Script to get the 1D histograms and save them as a PDF
    void Save1DHists(const char* print_name, const char* histname, std::string cut_type, bool scale);
    // -------------------------------------------------------------------------
    // Script to get the 2D histograms and save them as a PDF
    void Save2DHists(const char* print_name, const char* histname, std::string cut_type, bool yex);
    // -------------------------------------------------------------------------
    // Make plot of interaction efficiency
    void MakeInteractionEfficiency(const char *print_name);
    // -------------------------------------------------------------------------
    // Plot a 2d histogram normalised by row or column
    void Save2DHistsNorm(const char *print_name, const char *histname, std::string cut_type, bool yex, std::string normtype);
    // -------------------------------------------------------------------------


    private:


}; // End Class Histogram Plotter

#endif
