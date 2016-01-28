USE=config tools

# generic makefile
include ../maketools/make.module

PlumedMain.o: PlumedMainMap.tmp PlumedMainEnum.tmp

PlumedMainMap.tmp: PlumedMain.cpp
	../maketools/makecmd map < PlumedMain.cpp > PlumedMainMap.tmp

PlumedMainEnum.tmp: PlumedMain.cpp
	../maketools/makecmd enum < PlumedMain.cpp > PlumedMainEnum.tmp
	
