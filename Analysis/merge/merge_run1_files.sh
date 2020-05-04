# Command to merge the run1 files

# input mc file
mc_file=$1

# output merged file
outfile=$2

hadd -f -T $outfile files/nuexsec_data_run1.root files/nuexsec_dirt_run1.root files/nuexsec_ext_run1.root $mc_file
