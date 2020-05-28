LIBDIR = ../lib

# Compilation Flags
CXX = g++
CPPFLAGS = -I$(LIBDIR) -Wall -Wextra
CXXFLAGS=$(shell root-config --cflags)

# Libraries to inlcude e.g. ROOT
LIBS=$(shell root-config --libs) -L$(LIBDIR)

# Debug Flags
DEBUGFLAGS= -O0 -g -D _DEBUG

# Output name of the executable
MAIN = nuexsec

# List the object files here
# OBJ_FILES = bin/main.o bin/utility.o bin/selection.o bin/selection_cuts.o bin/Passed_Container.o bin/histogram_helper.o bin/histogram_plotter.o bin/SliceContainer.o bin/TreeHelper.o
OBJ_FILES = bin/*.o

# Main Target
all: $(MAIN)

$(MAIN): $(OBJ_FILES)
	@echo
	@echo -e  "\033[0;32mBuilding $(MAIN)\033[0;36m"
	$(CXX) $(LIBS) -o $@ $^
	@echo -e "\033[0;31mDone!"
	@mv include/*.gch bin
	@echo

# Build the Object files
bin/%.o: src/%.cxx include/%.h
	@mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEBUGFLAGS) -c $^
	@mv *.o bin

# Clean compilation
clean:
	@echo -e "\033[0;31m"
	rm nuexsec bin/* include/*.gch
	@echo
