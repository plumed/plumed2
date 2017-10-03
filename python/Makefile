# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean install

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the cython interface

else

ifdef python_bin

all:
	@echo Building python interface for PLUMED 
	CC="$(CC)" \
          LDSHARED="$(LDSHARED)" \
          CXX="$(CXX)" \
          program_name=plumed \
          plumedexe=../src/lib/plumed \
          include_dir=../src/wrapper \
          lib_dir=../src/lib \
          $(python_bin) buildPythonInterface.py build_ext -i

else

all:
	@echo Python support was not configure
endif

endif

clean:
	rm -fr *.so plumed.cpp build
