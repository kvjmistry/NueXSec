# Directory for TPC Object Libraries
LIBDIR = ../lib

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
OBJ_FILES = bin/main.o bin/utility.o bin/selection.o

# Main Target
all: $(MAIN)

$(MAIN): $(OBJ_FILES)
	echo Building $(MAIN)
	$(CXX) $(LIBS) -o $@ $^

# Build the Object files
bin/%.o: src/%.cxx include/%.h
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEBUGFLAGS) -c $^
	mv *.o bin

# Clean compilation
clean:
	rm nuexsec bin/*
