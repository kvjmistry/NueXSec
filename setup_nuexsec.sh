# Setup scrip for module to call necessary requirements of area

echo "Now compiling analysis script..."

echo "cd Analysis"
cd Analysis

echo "make"
make

echo "cd .."
cd ..

echo "setup complete!"
