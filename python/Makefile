# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

VERSION = $(shell if test -f ../VERSION ; then grep -v \\\# ../VERSION ; else echo "Unknown" ; fi )

.PHONY: all clean install

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the cython interface

else

ifdef python_bin

all:
	@echo Building python interface for PLUMED 
	unset CXX && unset CC && unset CFLAGS && unset CXXFLAGS && unset LDSHARED && \
          plumed_program_name=plumed \
          plumed_version=$(VERSION) \
          plumed_include_dir=../src/wrapper \
          plumed_default_kernel="$$PWD/../src/lib/libplumedKernel.$(SOEXT)" \
          $(python_bin) setup.py build_ext -i

else

all:
	@echo Python support was not configure
endif

endif

pip:
	cp ../src/wrapper/Plumed.h .
	cp ../VERSION .

clean:
	rm -fr *.so plumed.cpp build Plumed.h VERSION
