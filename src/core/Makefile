USE=config tools lepton small_vector mpiwrap

# generic makefile
include ../maketools/make.module

# -Wshadows helps to find slips in the code, 
# -Werror makes fixing warnings mandatory
CXXFLAGS:= $(CXXFLAGS) -Wshadow -Werror

CLEANLIST:=$(CLEANLIST) PlumedMainMap.inc PlumedMainEnum.inc GREXMap.inc GREXEnum.inc CLToolMainMap.inc CLToolMainEnum.inc ModuleMap.inc ActionDescriptions.inc

PlumedMain.o: PlumedMainMap.inc PlumedMainEnum.inc

PlumedMainMap.inc: PlumedMain.cpp
	../maketools/makecmd map < PlumedMain.cpp > PlumedMainMap.inc

PlumedMainEnum.inc: PlumedMain.cpp
	../maketools/makecmd enum < PlumedMain.cpp > PlumedMainEnum.inc

GREX.o: GREXMap.inc GREXEnum.inc

GREXMap.inc: GREX.cpp
	../maketools/makecmd map < GREX.cpp > GREXMap.inc

GREXEnum.inc: GREX.cpp
	../maketools/makecmd enum < GREX.cpp > GREXEnum.inc
	
CLToolMain.o: CLToolMainMap.inc CLToolMainEnum.inc

CLToolMainMap.inc: CLToolMain.cpp
	../maketools/makecmd map < CLToolMain.cpp > CLToolMainMap.inc

CLToolMainEnum.inc: CLToolMain.cpp
	../maketools/makecmd enum < CLToolMain.cpp > CLToolMainEnum.inc

ModuleMap.o: ModuleMap.inc

ModuleMap.inc: ../*/*.cpp $(wildcard ../../plugins/*/*.cpp ../../plugins/*/src/*.cpp)
	../maketools/makeModuleMap $@

ActionDescriptions.o: ActionDescriptions.inc

# FORCE makes the script rescan the directory so additions and removals are
# noticed as well as edits. The generator preserves the include's mtime when
# the extracted summaries have not changed.
ActionDescriptions.inc: FORCE $(wildcard ../*/*.cpp) ../maketools/makeActionDescriptions
	../maketools/makeActionDescriptions $@ ../..

.PHONY: FORCE
FORCE:
