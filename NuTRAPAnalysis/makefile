CXX := g++
RCINT := rootcint
RC := root-config

BDIR :=bin
LDIR :=lib

TOBJS := TransversityVariableObjects.cxx TransversityUtils.cxx GeneratorSpecifics.cxx
TOBJH := $(TOBJS:.cxx=.hxx)
TOBJDICTS := TransversityVariableObjects_dict.cxx
TDICTHEADERS := TransversityVariableObjects_dict.h
TOBJO := $(TOBJS:.cxx=.o)

UTILSLOC := ../utils
LIBUTILS := $(UTILSLOC)/lib/libPureGenUtils.a
LIBUTILSLD := -L$(UTILSLOC)/lib -lPureGenUtils

TARGET := NuTRAPAnalysis.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs`

CXXFLAGS := -fPIC $(ROOTCFLAGS) -I$(UTILSLOC) -g -std=c++11
LDFLAGS := $(ROOTLDFLAGS) $(LIBUTILSLD)

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
	$(CXX) $(CXXFLAGS) $(TARGETSRC) -o $@ $(TOBJO) $(LDFLAGS)

GeneratorSpecifics.o: GeneratorSpecifics.cxx GeneratorSpecifics.hxx TransversityUtils.hxx TransversityVariableObjects.hxx
	$(CXX) $(CXXFLAGS) -c GeneratorSpecifics.cxx -o $@

TransversityVariableObjects.o: TransversityVariableObjects.cxx TransversityVariableObjects.hxx TransversityUtils.hxx
	$(CXX) $(CXXFLAGS) -c TransversityVariableObjects.cxx -o $@

TransversityUtils.o: TransversityUtils.cxx TransversityUtils.hxx
	$(CXX) $(CXXFLAGS) -c TransversityUtils.cxx -o $@

EventLoop.o: tools/EventLoop.cxx tools/EventLoop.hxx
	$(CXX) $(CXXFLAGS) $< -I. -c -o $@

$(BDIR)/EvDump.exe: tools/EvDump.cxx TransversityVariableObjects.cxx TransversityVariableObjects.hxx EventLoop.o $(LIBUTILS)
	$(CXX) $(CXXFLAGS) tools/EvDump.cxx -I. -Itools -o EvDump.exe EventLoop.o $(TOBJO) $(LDFLAGS)
	mv EvDump.exe $(BDIR)

$(BDIR)/MomCons.exe: tools/MomCons.cxx TransversityVariableObjects.cxx TransversityVariableObjects.hxx EventLoop.o $(LIBUTILS)
	$(CXX) $(CXXFLAGS) tools/MomCons.cxx -I. -Itools -o MomCons.exe EventLoop.o $(TOBJO) $(LDFLAGS)
	mv MomCons.exe $(BDIR)

$(BDIR)/TauCalc.exe: tools/TauCalc.cxx TransversityVariableObjects.cxx TransversityVariableObjects.hxx EventLoop.o $(LIBUTILS)
	$(CXX) $(CXXFLAGS) tools/TauCalc.cxx -I. -Itools -o TauCalc.exe EventLoop.o $(TOBJO) $(LDFLAGS)
	mv TauCalc.exe $(BDIR)

$(BDIR)/EInvisCalc.exe: tools/EInvisCalc.cxx TransversityVariableObjects.cxx TransversityVariableObjects.hxx EventLoop.o $(LIBUTILS)
	$(CXX) $(CXXFLAGS) tools/EInvisCalc.cxx -I. -Itools -o EInvisCalc.exe EventLoop.o $(TOBJO) $(LDFLAGS)
	mv EInvisCalc.exe $(BDIR)

tools: $(BDIR)/TauCalc.exe $(BDIR)/EInvisCalc.exe


clean:
	rm -f $(TOBJDICTS)\
        $(TDICTHEADERS)\
        $(TOBJO)\
        $(TARGET)

distclean: clean
	rm -rf $(BDIR)
	cd $(UTILSLOC); $(MAKE) clean
