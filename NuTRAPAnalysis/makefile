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

UTILSLOC=../utils
LIBUTILS=libPureGenUtils.so
UTILSDEPSO:=$(UTILSLOC)/lib/$(LIBUTILS)

TARGET := NuTRAPAnalysis.exe
TARGETSRC := $(TARGET:.exe=.cxx)

ROOTCFLAGS := `$(RC) --cflags`
ROOTLDFLAGS := `$(RC) --libs --glibs`

CXXFLAGS := -fPIC $(ROOTCFLAGS) -g -std=c++11
LDFLAGS := $(ROOTLDFLAGS) -Wl,-rpath=.

.PHONY: all wtests clean

all: $(TARGET)
	mkdir -p $(BDIR) $(LDIR)
	mv $(TARGET) $(BDIR)/
	mv $(LIBUTILS) $(LDIR)/
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built NuTRAPAnalysis."
	@echo "*********************************************************************"

wtests: $(TARGET) tests
	./TransversityTests.exe
	mkdir -p $(BDIR) $(LDIR)
	mv $(TARGET) $(BDIR)/
	mv $(LIBUTILS) $(LDIR)/
	@echo ""
	@echo "*********************************************************************"
	@echo "Success. Built NuTRAPAnalysis."
	@echo "*********************************************************************"

$(LIBUTILS):
	cd $(UTILSLOC); $(MAKE)
	ln -s `readlink -f $(UTILSDEPSO)` $@

tests: TransversityUnitTests.cxx TransversityUtils.o
	$(CXX) -o TransversityTests.exe $(CXXFLAGS) $(LDFLAGS) TransversityUnitTests.cxx TransversityUtils.o

$(TARGET): $(TARGETSRC) $(TOBJO) $(LIBUTILS)
	$(CXX) -o $@ -I$(UTILSLOC) $(CXXFLAGS) $(LDFLAGS) $(TARGETSRC) $(TOBJO) $(LIBUTILS)

TransversityVariableObjects_dict.o: TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h TransversityUtils.o
	$(RCINT) -f TransversityVariableObjects_dict.cxx -c -p TransversityVariableObjects.hxx TransversityVariableObjects_linkdef.h
	$(CXX) $(CXXFLAGS) -c TransversityVariableObjects_dict.cxx -o $@

TransversityVariableObjects.o: TransversityVariableObjects.cxx TransversityVariableObjects.hxx TransversityUtils.o
	$(CXX) $(CXXFLAGS) -c TransversityVariableObjects.cxx -o $@

TransversityUtils.o: TransversityUtils.cxx TransversityUtils.hxx
	$(CXX) $(CXXFLAGS) -c TransversityUtils.cxx -o $@

clean:
	rm -rf $(TOBJDICTS)\
        $(TDICTHEADERS)\
        $(TOBJO)\
        $(LIBUTILS)\
        $(TARGET)\
        $(BDIR) \
        $(LDIR)
	cd $(UTILSLOC); $(MAKE) clean
