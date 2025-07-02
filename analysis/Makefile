#
#

CXXFLAGS = -g -O2 -Wall -Wuninitialized -fPIC

# required ZLIB library

CXXFLAGS += -DHAVE_LIBZ



# ROOT libraries

ifdef ROOTSYS

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTFEATURES := $(shell root-config --features)
ROOTGLIBS    := $(shell root-config --glibs) -lThread -lSpectrum
HAVE_ROOT_HTTP := $(findstring http,$(ROOTFEATURES))
HAVE_ROOT_XML  := $(findstring xml,$(ROOTFEATURES))

CXXFLAGS  += -DHAVE_ROOT $(ROOTCFLAGS)

ifdef HAVE_ROOT_XML
CXXFLAGS  += -DHAVE_ROOT_XML
ROOTGLIBS += -lXMLParser
endif

ifdef HAVE_ROOT_HTTP
CXXFLAGS  += -DHAVE_ROOT_HTTP -DHAVE_THTTP_SERVER
ROOTGLIBS += -lRHTTP
endif

endif # ROOTSYS


all: ANAISplusLevel2.exe ANAISplusLevel1.exe 

ANAISplusLevel2.exe: ANAISplusLevel2.cxx 
	$(CXX) -o $@ $(CXXFLAGS) $^ $(ROOTGLIBS) -lm -lz -lpthread -lutil

ANAISplusLevel1.exe: ANAISplusLevel1.cxx 
	$(CXX) -o $@ $(CXXFLAGS)  $^ $(ROOTGLIBS) -lm -lz -lpthread -lutil

%.o: %.cxx
	$(CXX) -o $@ $(CXXFLAGS) $(INCFLAGS) -c $<

dox:
	doxygen

clean::
	-rm -f *.o *.a *.so
	-rm -f *.exe
	-rm -f *Dict.cxx
	-rm -rf *.pcm

# end
