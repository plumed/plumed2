# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../Makefile.conf
endif

VERSION = $(shell if test -f ../VERSION.txt ; then grep -v \\\# ../VERSION.txt ; else echo "Unknown" ; fi )

.PHONY: all clean pypi pip

plumed_compiled := $(wildcard ../src/lib/plumed)

ifdef plumed_found_python_build

# The interface is rebuilt only if something it depends on has changed.
# It does not depend on the content of the kernel library, but only on its path
# (which is hardcoded as default kernel), so it does not need to be rebuilt when
# the C++ code changes.
# The stamp file records when the last build was completed.
# The config file records the settings used for the build (location of the tree,
# python executable, extension of shared libraries and version) and is only
# updated when they change.
STAMP := .python-build.stamp
CONFIG := .python-build.config

all: $(STAMP)

.PHONY: FORCE
FORCE:

$(CONFIG): FORCE
	@echo "$(CURDIR) $(python_bin) $(SOEXT) $(VERSION)" > $@~
	@cmp -s $@~ $@ || cp $@~ $@
	@rm -f $@~

$(STAMP): plumed.c plumed.pyx cplumed.pxd setup.py onlyBuild.sh README.rst MANIFEST.in \
          ../src/wrapper/Plumed.h ../Makefile.conf PLUMED_VERSION $(CONFIG)
	@echo Building python interface for PLUMED 
	rm -fr dist build *.so plumed.egg-info
	unset CXX && unset CC && unset CFLAGS && unset CXXFLAGS && unset LDSHARED && \
	  ./onlyBuild.sh -v $(VERSION) -p $(python_bin) \
          -k "$$PWD/../src/lib/libplumedKernel.$(SOEXT)"
	touch $@

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


pypi: clean PLUMED_VERSION plumed.c
	@echo "Setup the source distribution with the settings for pypi"
	$(python_bin) -m build --sdist .
	@echo "now use: $(python_bin) -m twine upload dist/plumed-$(VERSION).tar.gz"

pip: PLUMED_VERSION

# plumed.c is generated from plumed.pyx and cplumed.pxd with cython and it is committed,
# so that building the interface (make, make install, pip install) does not require cython
# and does not modify any tracked file.
# If you modify plumed.pyx or cplumed.pxd, plumed.c is regenerated here and
# must be committed as well. Use the reference cython version
# (it is checked in the CI that plumed.c is consistent with plumed.pyx).
# If cython is not available, the existing plumed.c is used, with a warning.
# "make cythonize" regenerates plumed.c unconditionally.
CYTHON_VERSION := 3.2.4

.PHONY: cythonize
cythonize:
	$(python_bin) -m cython plumed.pyx

plumed.c: cplumed.pxd plumed.pyx
	@if $(python_bin) -c "import Cython" 2>/dev/null ; then \
	  v="$$($(python_bin) -c 'import Cython; print(Cython.__version__)')" ; \
	  echo "Regenerating plumed.c with cython $$v" ; \
	  test "$$v" = "$(CYTHON_VERSION)" || echo "WARNING: the reference cython version is $(CYTHON_VERSION), plumed.c will differ from the one generated in the CI" ; \
	  $(python_bin) -m cython plumed.pyx ; \
	else \
	  echo "WARNING: plumed.pyx or cplumed.pxd are newer than plumed.c, but cython is not available with $(python_bin)" ; \
	  echo "WARNING: using the existing plumed.c, which might be out of date" ; \
	  echo "WARNING: install cython $(CYTHON_VERSION) and run make -C python cythonize to regenerate it" ; \
	  touch plumed.c ; \
	fi

clean:
	rm -fr dist build *.so plumed.egg-info PLUMED_VERSION *.log .python-build.stamp .python-build.config
