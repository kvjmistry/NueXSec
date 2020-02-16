# NueXSec
A repository for a Nue XSec Selection using the NuMI Beam

## Input to Code
This repository uses the input n-tuples written by the pandora LEE team. Please see this link for more info. I plan on integrating a numi branch to this code to allow for numi specific items. 

https://github.com/ubneutrinos/searchingfornues

## Analysis and Test Scripts
Be sure to read the comments in the relevant files you're using, the header files contain information what each class and function is doing!

MCC9:  

The analysis code lives in the Analysis folder. To compile, type `make` in the directory. 

To run you can do:
```
./nuexsec --mc path_to_mc_file --ext path_to_ext_file -- data path_to_on_beam_data_file
```

Full options are available by doing `./nuexsec --h`

The mc file is really an overlay file. I haven't built the functionality for fullmc yet.

Additional functionality includes providing the path to a custom configuration text file using `-c config.txt` (see utility.h for the formatting of the `configure` function used to set the cut values). If this is not set, the code will use a set of default parameters in `main.h`.

Also, the default running condition is to use the full selection and produce many many plots - use `--slim` to run more quickly and simply see the cut's performances.


## Calculating the Cross Section

To Be Updated....


## Other Modules

FlashValidation contains a LArSoft Module for analysing the flash information in the events. EventRate contains a LArSoft module to make the eventrate distribution. If you need to use these then ask me and I will document this further. 
