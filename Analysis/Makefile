# Directory for TPC Object Libraries
LIBDIR = ../lib
MODULESDIR = ../Modules

# Compilation Flags
CXX = g++
CPPFLAGS = -I$(LIBDIR) -Wall
CXXFLAGS=$(shell root-config --cflags)

# Libraries to inlcude e.g. ROOT, TPC Object Maker
LIBS=$(shell root-config --libs) -L$(LIBDIR) -lNueXSec

# Debug Flags
DEBUGFLAGS= -O0 -g -D _DEBUG

# Output name of the executable
MAIN = nuexsec

# List the object files here
OBJ_FILES = bin/main.o bin/utility.o bin/selection.o bin/selection_cuts.o bin/Passed_Container.o

# Main Target
all: $(MAIN)

$(MAIN): $(OBJ_FILES)
	@echo
	@echo -e  "\033[0;32mBuilding $(MAIN)\033[0;36m"
	$(CXX) $(LIBS) -o $@ $^
	@echo -e "\033[0;31mDone!"
	@echo

# Build the Object files
bin/%.o: src/%.cxx
	@mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEBUGFLAGS) -c $^
	@mv *.o bin

# Clean compilation
clean:
	@echo -e "\033[0;31m"
	rm nuexsec bin/*
	@echo
