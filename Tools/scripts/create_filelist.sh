
# First create a file list
source filelist_from_runsubrunlist.sh numi_uboone_run1_beamon_offset1_mcc9_reco2_v08_00_00_28_beam_good_bad ../../Analysis/files/run1_run_subrun_list_data.txt

# next create a sam def from the files for prestageing
source filelist_to_samdef.sh  kmistry_testdef_april2020 ../bin/files_run1_run_subrun_list_data.txt 

cd ../bin

# now create a filelist with the full path
samdef_to_filelist kmistry_testdef_april2020

cd ../scripts
