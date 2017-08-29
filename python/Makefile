# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../../Makefile.conf
endif

.PHONY: all clean

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the cython interface

else

all:
	python buildPythonInterface.py build_ext -i

endif

clean:
	rm -fr *.so plumed.cpp build
