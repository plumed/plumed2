# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

VERSION = $(shell if test -f ../VERSION ; then grep -v \\\# ../VERSION ; else echo "Unknown" ; fi )

.PHONY: all clean install

plumed_compiled := $(wildcard ../src/lib/plumed)

ifdef python_bin

all:
	@echo Building python interface for PLUMED 
	unset CXX && unset CC && unset CFLAGS && unset CXXFLAGS && unset LDSHARED && \
          plumed_program_name=plumed \
          plumed_force_cython=yes \
          plumed_version=$(VERSION) \
          plumed_include_dir=../src/wrapper \
          plumed_default_kernel="$$PWD/../src/lib/libplumedKernel.$(SOEXT)" \
          $(python_bin) setup.py build_ext -i

else

all:
	@echo Python support was not configure

endif

pip:
	cp ../VERSION .
	mkdir -p include
	cp ../src/wrapper/Plumed.h include/

ifndef python_bin
python_bin=python
endif

.PHONY: pypi

pypi:
	mkdir -p pypi
	mkdir -p pypi/include
	cp -r test pypi/
	cp ../src/wrapper/Plumed.h pypi/include/
	cp README.rst MANIFEST.in cplumed.pxd plumed.pyx setup.py ../VERSION pypi/
	cd pypi ; $(python_bin) setup.py sdist
	echo "now use: $(python_bin) -m twine upload dist/plumed-$(VERSION).tar.gz"
