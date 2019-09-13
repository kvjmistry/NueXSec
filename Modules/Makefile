CXX=g++
CXXFLAGS=$(shell root-config --cflags)
LIBS=$(shell root-config --libs)
LD=g++
LDFLAGS= -shared
DEBUGFLAGS= -O0 -D _DEBUG
COMPILE_FLAGS+=-fPIC

SOURCENAMES=TpcObjectContainer ParticleContainer
HEADERS=$(SOURCENAMES:%=%.h)
SOURCES=$(SOURCENAMES:%=%.cxx)
OBJECTS=$(SOURCENAMES:%=%.o)

DICTSOURCE=$(SOURCES)
TARGETNAME=NueXSec

DICTNAME=$(TARGETNAME)_dict
DICTOBJECTS=$(DICTNAME).o $(DICTNAME).cxx $(DICTNAME)_rdict.pcm
TARGETS=lib$(TARGETNAME).so lib$(TARGETNAME)_dict.so


all: $(TARGETS)

lib$(TARGETNAME)_dict.so: lib$(TARGETNAME).so
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

lib$(TARGETNAME).so: $(DICTNAME).cxx $(SOURCES) 
	$(LD) $(LDFLAGS) $(COMPILE_FLAGS) $(LIBS) $(CXXFLAGS) -o $@ $^


$(DICTNAME).cxx: $(HEADERS) LinkDef.h
	rootcling -f $@ -c $(CXXFLAGS) -p $^

debug:
	@echo "Target name: '$(TARGETNAME)'"

clean:
	rm -f $(OBJECTS) $(TARGETS) $(DICTOBJECTS)

copy:
	mv *NueXSec* ../lib/
