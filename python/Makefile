# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

VERSION = $(shell if test -f ../VERSION.txt ; then grep -v \\\# ../VERSION.txt ; else echo "Unknown" ; fi )

.PHONY: all clean pypi cythonize pip

plumed_compiled := $(wildcard ../src/lib/plumed)

ifdef plumed_found_python_build

all: clean PLUMED_VERSION
	@echo Building python interface for PLUMED 
	unset CXX && unset CC && unset CFLAGS && unset CXXFLAGS && unset LDSHARED && \
	  ./onlyBuild.sh -v $(VERSION) -p $(python_bin) \
          -k "$$PWD/../src/lib/libplumedKernel.$(SOEXT)" \

else

all:
	@echo Python support was not configured

endif
PLUMED_VERSION: ../VERSION.txt
	cp $< $@
#or ?grep -v '\#' $< $@

ifndef python_bin
python_bin=python
endif


pypi: clean PLUMED_VERSION
	@echo "Setup the source distribution with the settings for pypi"
	$(python_bin) -m build --sdist .
	@echo "now use: $(python_bin) -m twine upload dist/plumed-$(VERSION).tar.gz"

pip: PLUMED_VERSION

cythonize: plumed.c

plumed.c: cplumed.pxd plumed.pyx
	$(python_bin) -m cython plumed.pyx

clean:
	rm -fr dist build *.so plumed.egg-info PLUMED_VERSION *.log
