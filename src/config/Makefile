
# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../../Makefile.conf
endif

# sed script to edit program name defaults to do nothing
ifndef program_transform_name
  program_transform_name=s,x,x,
endif

.PHONY: all obj
ifdef PLUMED_FORCE
.PHONY: ConfigInstall.inc Config.inc version.h
endif

# default target
all:
	$(MAKE) PLUMED_FORCE=yes ConfigInstall.inc Config.inc version.h config.txt ModuleBuild.inc
	$(MAKE) Config.o ConfigInstall.o compile_options.sh

obj: all

config.txt: ../../Makefile.conf
	@echo "Updating config.txt file"
	@makefile_conf="$(realpath ../../Makefile.conf)" \
         python_bin="$(python_bin)" \
         mpiexec="$(mpiexec)" \
         ../maketools/update-config-txt.sh config.txt $(CPPFLAGS)

version.h:
	@echo "Updating version.h"
	@ ../maketools/update-version.sh version.h

ModuleBuild.inc: FORCE
	../maketools/module-build-info

.PHONY: FORCE
FORCE:

ConfigInstall.inc: Config.inc.in
	@echo "Updating ConfigInstall.inc"
	@sed "s/@SOEXT@/$(SOEXT)/g" Config.inc.in | \
         sed "s/@ISINSTALLED@/true/g" | \
         sed "s|@PLUMED_ROOT@|$(libdir)/$(program_name)|g" | \
         sed "s|@htmldir@|$(htmldir)|g" | \
         sed "s|@includedir@|$(includedir)|g" | \
         sed "s|@program_name@|$(program_name)|g" | \
         sed "s|@libdir@|$(libdir)|g" > $@~
	cmp -s $@~ $@ || cp $@~ $@
	rm $@~

# file to import compilation options inside c++
Config.inc: Config.inc.in
	@echo "Updating Config.inc"
	@sed "s/@SOEXT@/$(SOEXT)/g" $< | \
         sed "s/@ISINSTALLED@/false/g" | \
         sed "s|@htmldir@|xxxxNAxxxx|g" | \
         sed "s|@includedir@|xxxxNAxxxx|g" | \
         sed "s|@program_name@|xxxxNAxxxx|g" | \
         sed "s|@PLUMED_ROOT@|$(CURDIR:src/config=)|g" > $@~
	cmp -s $@~ $@ || cp $@~ $@
	rm $@~

links: 
	@echo Nothing to do in $(CURDIR)


Config.o: ModuleBuild.inc Config.inc Makefile.conf.xxd version.h
ConfigInstall.o: ModuleBuild.inc ConfigInstall.inc Makefile.conf.xxd version.h

Makefile.conf.xxd: ../../Makefile.conf
	cat ../../Makefile.conf | ../maketools/xxd > Makefile.conf.xxd

# file to import compilation options inside a bash script
compile_options.sh: ../../Makefile.conf Makefile
	@echo Building compile_options.sh, a sourceable bash script with some compilation options
	@echo "mklib_cxx='"$(CXX)"'" > $@
	@echo "compile='"$(CXX) -c -I'"$$PLUMED_INCLUDEDIR/$$PLUMED_PROGRAM_NAME"'/ $(CPPFLAGS) $(CXXFLAGS)"'" >> $@
# on osx, in order to use two-level namespace, we should back link the Kernel library
# the path to the kernel is different for installed and not-installed plumed.
# notice that $(DYNAMIC_LIBS) are needed as well in order to allow using functions
# that are exported from plumed (e.g. blas/lapack/mpi functions).
# on linux, the kernel is needed when it is loaded with RTLD_LOCAL (as done by
# the python interface), since its symbols are then not visible to the libraries
# loaded later. When the kernel is loaded globally (e.g. by the plumed executable)
# it is not necessary, due to flat namespace.
# As an alternative, on osx we could use -undefined dynamic_lookup. I prefer to use
# the default option (two level namespace) since it should give linking errors in advance,
# but we might change it later.
# mklib uses link_shared, followed by the objects and then by link_libs_(un)installed:
# the libraries must follow the objects using them, otherwise linkers using
# --as-needed (the default on some distributions, e.g. Ubuntu) drop them.
	@echo "link_shared='"$(LDSHARED)"'" >> $@
	@echo "link_libs_installed='"$(DYNAMIC_LIBS) $(libdir)/lib$(program_name)Kernel.$(SOEXT)"'" >> $@
	@echo "mklib_link_libs_installed='"$(filter -fopenmp -qopenmp -pthread -stdlib=%,$(CXXFLAGS)) -L$(libdir) -Wl,-rpath,$(libdir) $(libdir)/lib$(program_name)Kernel.$(SOEXT)"'" >> $@
# this is removed with a grep statement in src/lib/Makefile (see https://github.com/plumed/plumed2/issues/408)
	@echo "link_libs_uninstalled='"$(DYNAMIC_LIBS) $(CURDIR:/config=)/lib/libplumedKernel.$(SOEXT)"'" >> $@
# kept for scripts written for previous versions (libraries before the objects, see above)
	@echo "link_installed='"$(LDSHARED) $(DYNAMIC_LIBS) $(libdir)/lib$(program_name)Kernel.$(SOEXT)"'" >> $@
# this is removed with a grep statement in src/lib/Makefile (see https://github.com/plumed/plumed2/issues/408)
	@echo "link_uninstalled='"$(LDSHARED) $(DYNAMIC_LIBS) $(CURDIR:/config=)/lib/libplumedKernel.$(SOEXT)"'" >> $@
	@echo "mpi_provider_compile_flags='"$(CFLAGS)"'" >> $@
	@echo "mpi_provider_link_flags='"$(filter-out $(CXX),$(LDSHARED)) $(LDFLAGS)"'" >> $@
	@echo "soext='"$(SOEXT)"'" >> $@

.PHONY: clean
clean:
	rm -f Config.inc compile_options.sh Makefile_conf version.h config.txt ModuleBuild.inc
	rm -fr deps 
	rm -f Config.inc ConfigInstall.inc compile_options.sh install.conf
	rm -f $(CLEANLIST)

# generic makefile rules
include ../maketools/make.rules
