echo "generating FHC files..."
root -l -q -b 'Weight_Histograms.C("fhc")'
echo
echo "generating RHC files..."
root -l -q -b 'Weight_Histograms.C("rhc")'