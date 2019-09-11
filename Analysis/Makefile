LIBDIR = ../lib
CXX = g++
CXXFLAGS=$(shell root-config --cflags) -fPIC
LIBS=$(shell root-config --libs) -L$(LIBDIR) -lNueXSec 
DEBUGFLAGS= -O0 -g -D _DEBUG 
MAIN = nuexsec
SRC = src/*.cxx

all: $(MAIN)

$(MAIN) : $(SRC)
	@echo Building $(MAIN)
	@$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) $(LIBS) -o $@ $^

clean:
	rm nuexec