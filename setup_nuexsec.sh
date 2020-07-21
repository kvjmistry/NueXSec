# Setup scrip for module to call necessary requirements of area

echo "Now compiling analysis script..."

echo "cd Analysis"
cd Analysis

echo "make"
make

echo "cd Systematics"
cd Systematics

echo "source generate_flux_weight_hists.sh"
source generate_flux_weight_hists.sh

echo "cd ../.."
cd ../..

echo "setup complete!"
