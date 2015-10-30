CXX := g++
RCINT := rootcint
RC := root-config

BDIR :=bin
LDIR :=lib

TOBJS := TransversityVariableObjects.cxx TransversityUtils.cxx
TOBJH := $(TOBJS:.cxx=.hxx)
TOBJDICTS := TransversityVariableObjects_dict.cxx
TDICTHEADERS := TransversityVariableObjects_dict.h
TOBJO := $(TOBJS:.cxx=.o)
TOBJO += TransversityVariableObjects_dict.o

UTILSLOC := ../utils
LIBUTILS := $(UTILSLOC)/lib/libPureGenUtils.a
LIBUTILSLD := -L$(UTILSLOC)/lib -lPureGenUtils

TARGET := NuTRAPAnalysis.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs`

CXXFLAGS := -fPIC $(ROOTCFLAGS) -g -std=c++11
LDFLAGS := $(ROOTLDFLAGS) $(LIBUTILSLD) -Wl,-rpath=.

.PHONY: all wtests clean

all: $(TARGET)
	mkdir -p $(BDIR)
	mv $(TARGET) $(BDIR)/
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built NuTRAPAnalysis."
	@echo "*********************************************************************"

wtests: $(TARGET) tests
	./TransversityTests.exe
	mkdir -p $(BDIR)
	mv $(TARGET) $(BDIR)/
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built NuTRAPAnalysis."
	@echo "*********************************************************************"

$(LIBUTILS):
	cd $(UTILSLOC); $(MAKE)

tests: TransversityUnitTests.cxx TransversityUtils.o
	$(CXX) -o TransversityTests.exe $(CXXFLAGS) TransversityUnitTests.cxx TransversityUtils.o $(LDFLAGS)

$(TARGET): $(TARGETSRC) $(TOBJO) $(LIBUTILS)
	$(CXX) $(CXXFLAGS) $(TARGETSRC) -o $@ -I$(UTILSLOC) $(TOBJO) $(LDFLAGS)

TransversityVariableObjects_dict.o: TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h TransversityUtils.hxx
	$(RCINT) -f TransversityVariableObjects_dict.cxx -c -p TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h
	$(CXX) $(CXXFLAGS) -c TransversityVariableObjects_dict.cxx -o $@

TransversityVariableObjects.o: TransversityVariableObjects.cxx TransversityVariableObjects.hxx TransversityUtils.hxx
	$(CXX) $(CXXFLAGS) -c TransversityVariableObjects.cxx -o $@

TransversityUtils.o: TransversityUtils.cxx TransversityUtils.hxx
	$(CXX) $(CXXFLAGS) -c TransversityUtils.cxx -o $@

clean:
	rm -f $(TOBJDICTS)\
        $(TDICTHEADERS)\
        $(TOBJO)\
        $(TARGET)

distclean: clean
	rm -rf $(BDIR) \
	cd $(UTILSLOC); $(MAKE) clean
