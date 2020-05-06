# Command to merge the run 3 files

# input mc file
mc_file=$1

# output merged file
outfile=$2

# Command to merge the run3 files
hadd -f -T $outfile files/nuexsec_data_run3.root files/nuexsec_dirt_run3.root files/nuexsec_ext_run3.root $mc_file
