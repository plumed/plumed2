# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

.PHONY: all clean install

plumed_compiled := $(wildcard ../src/lib/plumed)

ifeq ($(strip $(plumed_compiled)),)

all:
	@echo You must compile plumed before building the cython interface

install:
	@echo You must compile plumed before installing the cython interface

else

ifneq (,$(findstring __PLUMED_HAS_PYTHON,$(CPPFLAGS)))

all:
	python buildPythonInterface.py build_ext -i

install:
	python buildPythonInterface.py install

else

all:
	@echo Did not find __PLUMED_HAS_PYTHON
	@echo $(CPPFLAGS)

install:

endif

endif

clean:
	rm -fr *.so plumed.cpp build
