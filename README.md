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
https://github.com/kvjmistry/searchingfornues/tree/numi

## POT Counting
To get the POT and number of trigers from the pandora ntuples, go to the POT directory and do `source run_GetPOT.sh <input file> <mc/data>`

if data file is given then zarko's script will be called into action.

Please see this wiki for how to read these numbers:
https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/NuMI_Documentation#How-to-do-NuMI-POT-counting

Make sure to update the POT numbers in the main.h header when running the selection next.


## Analysis and Test Scripts
Be sure to read the comments in the relevant files you're using, the header files contain information what each class and function is doing!

MCC9:  

The analysis code lives in the Analysis folder. To compile, type `make` in the directory. 

To run you can do:
```
./nuexsec --run <1/3b> [--mc path_to_overlay_file] [--ext path_to_ext_file] [--data path_to_on_beam_data_file] [--dirt path_to_dirt_overlay_file]
```

Full options are available by doing `./nuexsec --h`

Note that mc is configuring an overlay file, same for dirt.

Configuration quantities such as the POT and cut values can be found in `main.h`. I will eventually move the POT config to a text file when there are detector variations. Still need to work on confiuring all the cuts from this file. 

Also, the default running condition is to use the full selection and produce many many plots - use `--slim` to run more quickly and simply see the cut's performances.

The output of the selection will be in the files directory. We make 1 file for each type (mc/dirt/ext/data). To make the full plots we combine these files with `hadd merged_file_run1.root *run1*`. We can then run the histogram maker:
`./nuexsec --run 1 --hist <merged_file_run1.root>`

I will try to make some scripts to make the merging more seamless once the selection is in less of a chaotic state.

## Other Modules

FlashValidation contains a LArSoft Module for analysing the flash information in the events. EventRate contains a LArSoft module to make the eventrate distribution. If you need to use these then ask me and I will document this further. 
