# NueXSec
A repository for a Nue XSec Selection using the NuMI Beam

## Setup instructions
Clone the repository by doing `git clone https://github.com/kvjmistry/NueXSec.git`
If you want the current branch I am working on, then you can do `git checkout pandora`

next source the setup script: `source setup_nuexsec.sh`

and thats it!

What follows is whats in this repository and how to use it.

## Input to Code
This repository uses the input n-tuples written by the pandora LEE team. Please see this link for more info. 

https://github.com/ubneutrinos/searchingfornues

I have forked this repository from the v30branch to make some numi chages. You can find the modified code in the numi branch of this fork:
https://github.com/kvjmistry/searchingfornues/tree/numi_tag_022820

## POT Counting
To get the POT and number of trigers from the pandora ntuples, go to the POT directory and do `source run_GetPOT.sh <input file> <mc/data>`

if data file is given then zarko's script will be called into action.

Please see this wiki for how to read these numbers:
https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/NuMI_Documentation#How-to-do-NuMI-POT-counting

Make sure to update the POT numbers in the main.h header when running the selection next.


## Running the Selection Code
Be sure to read the comments in the relevant files you're using, the header files contain information what each class and function is doing!

Compilation of the code can be done by simply entering `make` in the Analysis folder. This will create a binary nuexsec which we can use to run the code. If you make a change to the code, simply recompile by entering the `make` command. Sometimes you can get an error, so recompiling from scratch can help. To do this do `make clean` followed by `make`.

The analysis code lives in the Analysis folder. To compile, type `make` in the directory. 

A full list of commands that are being run in the analysis is listed in the run_selection_mcc9_run1.sh files. There are for running on Krish's local computer. If you want to run on the gpvm, then there should be a `_gpvm` version of the file.

In this instructions, I show you how to run over Run1. Hopefully its obvious how you can adapt this to run Run3.

To run the selection over a searchingfornues NuMI ntuple, you can do:
```
./nuexsec --run 1 [--mc path_to_overlay_file] 
./nuexsec --run 1 [--ext path_to_ext_file]
./nuexsec --run 1 [--data path_to_on_beam_data_file]
./nuexsec --run 1 [--dirt path_to_dirt_overlay_file]
```
This will produce a histogram file in the `files` folder with the output of the selection at each selection stage (amongst other plots) and another root file in `files/trees/` with variables that passed the selection. One file for mc, dirt, ext and data is produced in each of these cases. 

Note that mc is configuring an overlay file, same for dirt.

In general, there are not as many nue events in the standard MC file. We can replace these events by running over the intrinsic nue file:

`./nuexsec --run 1 --mc [--mc path_to_intrinsic_nue_file]  --intrinsic intrinsic` 

# Printing the Selection
To print the selection results (which read the files in the `files/trees` folder) we can run the command:

`./nuexsec --run 1 --printonly --printall`
additional options for only printing data, ext etc are available too.
# Configuration and Options
Configuration quantities such as the POT and fiducial volume boundary can be found in `config.txt`. The Utility class `Initialise` function reads in these values along with other options you give. These parameters are then inherited across all classes in the code which you can use. 

Full options and more documentation is available by doing `./nuexsec --h`. 

# Making Histograms
To make plots from the output of the selection, we need to merge the individual files we produced in the selection stage. I have wrote a merge script to combine these files to one and can be run like:

`source merge/merge_run1_files.sh files/nuexsec_mc_run1.root files/nuexsec_run1_merged.root`

This merged file is what we use for the histogram plotter.

# Running the Cross Section Code

To run the cross section code, we need to run over the events that we have selected. Before doing this, we need to merge the TTrees into one file to make this easier:

`root -l -b -q 'merge/merge_uneaventrees.C("1", true, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "")'`

The bool is currently set to true. This is where we replace the nue events in the standard tree with the intrinsic ones. Set this to false to use the standard nue events (not as many stats!). 

Now we have the merged file, we can calculate the cross section by running:

`./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode default`

This will print a bunch of stuff and also make a file: `files/crosssec_run1.root` with the final output histograms.

# Running the Systematics and Reweighting Events

To run the systematics code, we run the cross section helper in `reweight` mode. We specify the type of systematic to run by and the variable to bin the cross section in:
```
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel unisim  --xsecvar elec_E
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel ppfx    --xsecvar elec_E
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel genie   --xsecvar elec_E
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel reint   --xsecvar elec_E
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel mcstats --xsecvar elec_E
```
This creates a new folder in `files/crosssec_run1.root` with the systematic plots you have just ran.

# Running the Systematics: Detector Variations
Detector variations are made in different samples so we cant use the CV file. We need to run the selection, create a unique file and then run the xsec code over it like it was a CV. E.g. in the case of `LY_Down`, this can be done by:

```
./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_LY_Down.root LY_Down

# Overwrite the true nue information
./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_LY_Down_intrinsic.root LY_Down --intrinsic intrinsic

# Merge and plot histograms
source merge/merge_run1_files.sh files/nuexsec_mc_run1_LY_Down.root files/nuexsec_run1_LY_Down_merged.root
./nuexsec --run 1 --hist files/nuexsec_run1_LY_Down_merged.root --var dummy LY_Down

root -l -b -q 'merge/merge_uneaventrees.C("1", true, "files/trees/nuexsec_selected_tree_mc_run1_'"LY_Down"'.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "'"LY_Down"'")'

# Run the xsec code
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_LY_Down.root --var dummy LY_Down --xsecmode default --xsecvar elec_E
```

Thankfully we have wrote an automation of this process which can be run by doing:
`source run_selection_mcc9_run1.sh var <variation name>`

or to run all of them in one go:
`source run_selection_mcc9_run1.sh allvar`

Note, this overwrites the CV name with `detvar_CV` in `files/crosssec_run1.root` so we can compare the variation to the right CV when we plot etc.

# Plotting the Systematics and Final Cross Section

We now need to plot the cross section and calculate the uncertainties along with the covariance matrices etc. You can do this by running 
`./nuexsec --run 1 --sys reweight --xsecvar elec_E`
This reads in all the stuff in `files/crosssec_run1.root`, and creates all the final plots.

# Creating Systematics by Selection Stage
This is pretty cpu intensive, so try to minise how often you run this. But basically to calculate the systematics for select plots at each cut stage, we can use a different mode of the cross section helper:
This is running the so called rw_cuts stage and adds a cut folder in `files/crosssec_run1.root`. *This can take hours to run*.
```
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel unisim --xsecplot rw_cuts --intrinsic intrinsic
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel ppfx   --xsecplot rw_cuts --intrinsic intrinsic
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel genie  --xsecplot rw_cuts --intrinsic intrinsic
./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel reint  --xsecplot rw_cuts --intrinsic intrinsic
```

# Creating File Lists for Events
We might want to save the selected events with their weights in a txt file so we can use with a package like ReMu or give to someone else to test. You can do this by running the cross section code in a different mode:
`./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode txtlist --xseclabel all --xsecvar elec_E`
This prints all systematics and weights to file (incl all universes!). If you want to just look at the unisim weights in each event, then use `--xseclabel unisim` etc.

## Other Modules

FlashValidation contains a LArSoft Module for analysing the flash information in the events. EventRate contains a LArSoft module to make the eventrate distribution. If you need to use these then ask me and I will document this further. 

## Tools

This folder contains scripts not directly related to the selection, but are useful. It currently contains a script for converting a run subrun list to a filelist which can be filtered on to get all the selected nue events.




