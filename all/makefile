CXX := g++
RCINT := rootcint
RC := root-config

BDIR :=bin
LDIR :=lib

TOBJS := TransversityVariableObjects.cxx TransversityUtils.cxx
TOBJH := $(TOBJS:.cxx=.hxx)
TOBJDICTS := $(TOBJS:.cxx=_dict.cxx)
TDICTHEADERS := $(TOBJS:.cxx=_dict.h)
TOBJLINKDEFS := $(TOBJS:.cxx=_linkdef.h)
TOBJSO := $(TOBJS:.cxx=.so)

UTILSLOC=../../utils
UTILSDEPSO=libPureGenUtils.so

TARGET := ProcessRooTrackerToTransversityVariables.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs`

CXXFLAGS := -fPIC $(ROOTCFLAGS) -g -std=c++11
LDFLAGS := $(ROOTLDFLAGS) -Wl,-rpath=.

.PHONY: all clean

all: $(TARGET) tests
	./TransversityTests.exe
	mkdir -p $(BDIR) $(LDIR)
	mv $(TARGET) $(BDIR)/
	mv $(UTILSDEPSO) $(TOBJSO) $(LDIR)/
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built ProcessRooTrackerToTransversityVariables."
	@echo "*********************************************************************"

$(UTILSDEPSO):
	ln -s `readlink -f $(UTILSLOC)/lib/$@` $@

tests: TransversityUnitTests.cxx $(TOBJH)
	$(CXX) -o TransversityTests.exe $(CXXFLAGS) $(LDFLAGS) TransversityUnitTests.cxx TransversityUtils.so

$(TARGET): $(TARGETSRC) $(TOBJSO) $(TOBJH) $(UTILSDEPSO)
	$(CXX) -o $@ -I$(UTILSLOC) $(CXXFLAGS) $(LDFLAGS) $(TARGETSRC) $(TOBJSO) $(UTILSDEPSO)

TransversityVariableObjects_dict.cxx: TransversityVariableObjects.hxx TransversityUtils.hxx
	$(RCINT) -f $@ -c -p TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h

TransversityVariableObjects.so: TransversityVariableObjects.cxx TransversityVariableObjects_dict.cxx TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h TransversityUtils.hxx
	$(CXX) $(CXXFLAGS) -shared TransversityVariableObjects.cxx TransversityVariableObjects_dict.cxx -o $@

TransversityUtils.so: TransversityUtils.hxx
	$(CXX) $(CXXFLAGS) TransversityUtils.cxx -shared -o $@

clean:
	rm -rf $(TOBJDICTS)\
        $(TDICTHEADERS)\
        $(TOBJSO)\
        $(UTILSDEPSO)\
        $(TARGET)\
        $(BDIR) \
        $(LDIR)
