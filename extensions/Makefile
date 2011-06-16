
-include ../Makefile.conf

SRC=$(wildcard *.cpp)

OBJ=$(SRC:.cpp=.o)
SOBJ=$(SRC:.cpp=.$(SOEXT))

.PHONY: all clean

ifdef GCCDEP
all: libplumedExtensions.$(SOEXT) $(SOBJ)
else
all:
	@echo No configuration available
	@echo First run ./configure.sh in the root directory
endif


libplumedExtensions.$(SOEXT): $(OBJ)
	$(LDSO) -o $@ $^

clean:
	rm -f makefile.dep
	rm -f *.o
ifdef SOEXT
	rm -f *.$(SOEXT)
endif

makefile.dep : $(SRC) *.h
ifdef GCCDEP
	for i in $(ALL_SRC) ; do $(GCCDEP) $(CPPFLAGS) -MM "$${i}"; done > $@
endif

%.$(SOEXT) : %.o
	$(LDSO) -o $@ $^

-include makefile.dep

