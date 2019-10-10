# This will run the selection over the run1 files from the MCC8 era
#./nuexsec --mc /uboone/data/users/kmistry/work/NueXSection_Outputs/files/filter_numi_cosmic_v4.root --data /uboone/data/users/chill2/filter_data_numi_v7_full.root --dirt /uboone/data/users/chill2/filter_numi_dirt_v2_full.root --ext /uboone/data/users/chill2/filter_ext_v7_full.root

# Runs it in the background and supresses the terminal output
./nuexsec --mc /uboone/data/users/kmistry/work/NueXSection_Outputs/files/filter_numi_cosmic_v4.root &>/dev/null &
sleep $[ ( $RANDOM % 10 )  + 1 ]s
./nuexsec --data /uboone/data/users/chill2/filter_data_numi_v7_full.root &>/dev/null &
sleep $[ ( $RANDOM % 10 )  + 1 ]s
./nuexsec --ext /uboone/data/users/chill2/filter_ext_v7_full.root &>/dev/null &
sleep $[ ( $RANDOM % 10 )  + 1 ]s
./nuexsec --dirt /uboone/data/users/chill2/filter_numi_dirt_v2_full.root &>/dev/null &
