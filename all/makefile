CXX := g++
RCINT := rootcint
RC := root-config

TOBJS := TransversityVariableObjects.cxx
TOBJH := $(TOBJS:.cxx=.hxx)
UTILHEADER := TransversityUtils.hxx
TOBJDICTS := $(TOBJS:.cxx=_dict.cxx)
TDICTHEADERS := $(TOBJS:.cxx=_dict.h)
TOBJLINKDEFS := $(TOBJS:.cxx=_linkdef.h)
TOBJSO := $(TOBJS:.cxx=.so)

TARGET := ProcessRooTrackerToTransversityVariables.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs`

CXXFLAGS := -fPIC $(ROOTCFLAGS) -g
LDFLAGS := $(ROOTLDFLAGS) -Wl,-rpath=.

.PHONY: all clean

all: $(TARGET) tests
	TransversityTests.exe
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Note you will need to run me from here or add this directory"
	@echo "to your LD_LIBRARY_PATH."
	@echo "*********************************************************************"

tests: TransversityUnitTests.cxx $(UTILHEADER)
	$(CXX) -o TransversityTests.exe $(CXXFLAGS) $(LDFLAGS) TransversityUnitTests.cxx

$(TARGET): $(TARGETSRC) $(TOBJSO) $(UTILHEADER)
	$(CXX) -o $@ -I$(NEUTCLASSDIR) $(CXXFLAGS) $(LDFLAGS) $^

TransversityVariableObjects_dict.cxx: TransversityVariableObjects.hxx

TransversityVariableObjects.so: TransversityVariableObjects.cxx TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h $(UTILHEADER)
	$(RCINT) -f TransversityVariableObjects_dict.cxx -c -p TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h
	$(CXX) $(CXXFLAGS) -shared TransversityVariableObjects.cxx TransversityVariableObjects_dict.cxx -o $@

clean:
	rm -f $(TOBJDICTS)\
              $(TDICTHEADERS)\
              $(TOBJSO)\
              $(TARGET)\
