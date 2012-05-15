
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
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f *.tmp
	rm -f *.so
	rm -f *.dylib

%.$(SOEXT) : %.o
	$(LDSO) -o $@ $^

# this compiler is used for dependencies generation
# if unset, use normal CXX
ifndef CXXDEP
CXXDEP=$(CXX)
endif

%.o: %.cpp
	$(CXXDEP) -c -MM -MF$*.d $(CPPFLAGS) $(CXXFLAGS) $*.cpp
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	 sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp
ifndef XLF
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $*.cpp -o $*.o
endif

# Individual dependencies

ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_SRC_DEP)
endif

