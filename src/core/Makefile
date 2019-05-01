USE=config tools lepton

# generic makefile
include ../maketools/make.module

CLEANLIST:=$(CLEANLIST) PlumedMainMap.inc PlumedMainEnum.inc GREXMap.inc GREXEnum.inc CLToolMainMap.inc CLToolMainEnum.inc

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

