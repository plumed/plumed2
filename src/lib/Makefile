# here we build the list of modules to be linked in plumedKernel
KERNEL_MODULES := $(addprefix ../,$(shell ../maketools/find_modules.sh))

# This adds modules which are built anyway and linked with some ad hoc procedure
DIRS := $(sort ../config ../wrapper ../main  $(KERNEL_MODULES))

# Targets to build links
DIRSLINKS:=$(addsuffix .links,$(DIRS))

# include the machine dependent configuration
-include ../../Makefile.conf

# include installation info
-include ../config/install.conf

# if machine dependent configuration has been found:
ifdef GCCDEP

# source for plumed wrapper
SRC_WRAPPER=../wrapper/PlumedStatic.cpp

# source for plumed dynamic wrapper
SRC_DYNAMIC_WRAPPER=../wrapper/Plumed.c

# source for plumed main program
# this should not be linked into libraries
SRC_MAIN=../main/main.cpp

# source for plumed kernel
# everything else
SRC_KERNEL := $(shell for dir in $(KERNEL_MODULES) ; do ls $$dir/*.cpp 2>/dev/null ; ls $$dir/*.c 2>/dev/null ; done) ../config/Config$(PLUMED_INSTALL).o

# corresponding objects:
OBJ_WRAPPER := ../wrapper/PlumedStatic.o
OBJ_DYNAMIC_WRAPPER := ../wrapper/Plumed.o
OBJ_MAIN := $(SRC_MAIN:.cpp=.o)
XXX := $(SRC_KERNEL:.cpp=.o)
OBJ_KERNEL := $(XXX:.c=.o)
ALL_OBJ := $(OBJ_KERNEL) $(OBJ_WRAPPER) $(OBJ_DYNAMIC_WRAPPER) $(OBJ_MAIN)

# plumed kernel
PLUMED_KERNEL := $(PLUMED_LIBDIR)libplumedKernel.$(SOEXT)

# plumed library (kernel + wrapper)
PLUMED_SHARED_OBJ := $(PLUMED_LIBDIR)libplumed.$(SOEXT)

# plumed scripts
PLUMED_SCRIPTS_SH := $(shell cd ../../scripts/ ; ls *.sh 2>/dev/null)
PLUMED_SCRIPTS := $(PLUMED_SCRIPTS_SH:.sh=)
PLUMED_SCRIPTS_EXE := $(addprefix $(PLUMED_LIBDIR)plumed-,$(PLUMED_SCRIPTS))

# plumed main exec
PLUMED_MAIN_STATIC := $(PLUMED_LIBDIR)plumed
PLUMED_MAIN_SHARED := $(PLUMED_LIBDIR)plumed-shared
PLUMED_MAIN_RUNTIME := $(PLUMED_LIBDIR)plumed-runtime

ifdef PLUMED_LIBDIR_NODEST
  PLUMED_ROOT :="$(PLUMED_LIBDIR_NODEST)../../"
else
  PLUMED_ROOT :="$(realpath ../../)"
endif


.PHONY: all links $(DIRS) $(DIRSLINKS) all-static Plumed.cmake Plumed.inc
ifdef SOEXT
all:
	$(MAKE) dirslinks
	$(MAKE) dirs
	$(MAKE) lib
else
all:
	$(MAKE) all-static
endif

all-static:
	$(MAKE) dirslinks
	$(MAKE) dirs
	$(MAKE) lib-static

dirs:
	@echo "*** Compiling all directories ***"
	$(MAKE) $(DIRS)


install:
	cd ../config ; $(MAKE)
# two steps to be sure that ../config/install.conf is up to date
	$(MAKE) install-do

install-do:
	$(MAKE) uninstall
	mkdir -p $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/bin
	mkdir -p $(DESTDIR)$(PLUMED_INSTALL_ROOT)/src/lib
	mkdir -p $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/include/$(PLUMED_PROGRAM_NAME)
# rebuild libraries and main using new target dir
	$(MAKE) PLUMED_LIBDIR=$(DESTDIR)$(PLUMED_INSTALL_ROOT)/src/lib/ PLUMED_LIBDIR_NODEST=$(PLUMED_INSTALL_ROOT)/src/lib/ PLUMED_INSTALL=Install all
# copy .h, .cpp .o compile_options.sh and documentation:
	cd ../../ ; tar cf - patches scripts src/*/*.h src/*/*.cpp src/*/*.o src/config/compile_options.sh *-doc/html *-doc/regtests *-doc/manual.pdf | tar xf - -C $(DESTDIR)$(PLUMED_INSTALL_ROOT)/
# also copy .h files into include/ dir
	cd ../../src ; tar cf - */*.h */*/*.h | tar xf - -C $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/include/$(PLUMED_PROGRAM_NAME)/
# copy inc files changing the path
	for file in Plumed.inc* Plumed.cmake* ; do sed "s| $(realpath ../..)| $(PLUMED_INSTALL_ROOT)|g" $$file > $(DESTDIR)$(PLUMED_INSTALL_ROOT)/src/lib/$$file ; done
# modulefile
	@sed "s|@_PREFIX_@|$(PLUMED_INSTALL_PREFIX)|" modulefile.in | sed "s|@_SOEXT_@|$(SOEXT)|" | sed "s|@_PROGNAME_@|$(PLUMED_PROGRAM_NAME)|" > $(DESTDIR)$(PLUMED_INSTALL_ROOT)/src/lib/modulefile
# install links
	cd $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/bin/ && ln -s ../lib/$(PLUMED_PROGRAM_NAME)/src/lib/plumed $(PLUMED_PROGRAM_NAME)
	cd $(DESTDIR)$(PLUMED_INSTALL_ROOT) && for file in $(PLUMED_SCRIPTS) ; do ln -s src/lib/plumed-$$file plumed-$$file ; done
	cd $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/bin/ && ln -s ../lib/$(PLUMED_PROGRAM_NAME)/src/lib/plumed-patch $(PLUMED_PROGRAM_NAME)-patch
ifdef SOEXT
	cd $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/lib/ && ln -s $(PLUMED_PROGRAM_NAME)/src/lib/libplumed.$(SOEXT) lib$(PLUMED_PROGRAM_NAME).$(SOEXT)
	cd $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/lib/ && ln -s $(PLUMED_PROGRAM_NAME)/src/lib/libplumedKernel.$(SOEXT) lib$(PLUMED_PROGRAM_NAME)Kernel.$(SOEXT)
endif
	chmod -R go+rX,go-w $(DESTDIR)$(PLUMED_INSTALL_ROOT)
	chmod -R go+rX,go-w $(DESTDIR)$(PLUMED_INSTALL_ROOT)/plumed-*
	chmod -R go+rX,go-w $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/include/$(PLUMED_PROGRAM_NAME)
	chmod -R go+rX,go-w $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/bin/$(PLUMED_PROGRAM_NAME)
	@echo
	@echo "*** PLUMED has been installed ***"
	@echo
ifdef DESTDIR
	@echo "WARNING: Staged install with DESTDIR = $(DESTDIR)"
endif
	@echo "Install prefix : $(PLUMED_INSTALL_PREFIX)"
	@echo "Full name      : $(PLUMED_PROGRAM_NAME)"
	@echo
	@echo "Setup your environment"
	@echo "- Ensure this is in your execution path : $(PLUMED_INSTALL_PREFIX)/bin"
	@echo "- Ensure this is in your include path   : $(PLUMED_INSTALL_PREFIX)/include"
ifdef SOEXT
	@echo "- Ensure this is in your library path   : $(PLUMED_INSTALL_PREFIX)/lib"
	@echo "For runtime binding:"
	@echo "- Set this environment variable         : PLUMED_KERNEL=$(PLUMED_INSTALL_PREFIX)/lib/lib$(PLUMED_PROGRAM_NAME)Kernel.$(SOEXT)"
endif
	@echo
	@echo "To create a tcl module that sets all the variables above, use this one as a starting point:"
	@echo "$(PLUMED_INSTALL_ROOT)/src/lib/modulefile"
	@echo
	@echo "To uninstall, remove the following files and directories:"
	@echo $(PLUMED_INSTALL_ROOT)
	@echo $(PLUMED_INSTALL_PREFIX)/include/$(PLUMED_PROGRAM_NAME)
	@echo $(PLUMED_INSTALL_PREFIX)/bin/$(PLUMED_PROGRAM_NAME)
	@echo $(PLUMED_INSTALL_PREFIX)/bin/$(PLUMED_PROGRAM_NAME)-patch
ifdef SOEXT
	@echo $(PLUMED_INSTALL_PREFIX)/lib/lib$(PLUMED_PROGRAM_NAME).$(SOEXT)
	@echo $(PLUMED_INSTALL_PREFIX)/lib/lib$(PLUMED_PROGRAM_NAME)Kernel.$(SOEXT)
endif
ifeq ($(program_can_run),no)
	@echo "WARNING: $(PLUMED_PROGRAM_NAME) executable will not run on this machine"
	@echo "WARNING: to patch an MD code use 'plumed-patch'"
	@echo "WARNING: a few other command line tools are avalable as $(PLUMED_INSTALL_ROOT)/plumed-*"
else ifeq ($(program_can_run_mpi),no)
	@echo "WARNING: $(PLUMED_PROGRAM_NAME) executable will not run on this machine"
	@echo "WARNING: unless you invoke it as '$(PLUMED_PROGRAM_NAME) --no-mpi'"
	@echo "WARNING: to patch an MD code now use '$(PLUMED_PROGRAM_NAME) --no-mpi patch'"
	@echo "WARNING: (notice that MPI will be available anyway in the patched code)"
	@echo "WARNING: all command line tools are available as '$(PLUMED_PROGRAM_NAME) --no-mpi name-of-the-tool'"
	@echo "WARNING: e.g. '$(PLUMED_PROGRAM_NAME) --no-mpi driver'"
endif


uninstall:
	rm -fr $(DESTDIR)$(PLUMED_INSTALL_ROOT)
	rm -fr $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/include/$(PLUMED_PROGRAM_NAME)
	rm -f $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/bin/$(PLUMED_PROGRAM_NAME)
	rm -f $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/bin/$(PLUMED_PROGRAM_NAME)-patch
	rm -f $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/lib/lib$(PLUMED_PROGRAM_NAME).$(SOEXT)
	rm -f $(DESTDIR)$(PLUMED_INSTALL_PREFIX)/lib/lib$(PLUMED_PROGRAM_NAME)Kernel.$(SOEXT)

$(DIRS):
	$(MAKE) -C $@ obj

.PHONY: lib lib-static lib-shared

lib: lib-static lib-shared

lib-static: $(PLUMED_MAIN_STATIC) $(OBJ_WRAPPER) Plumed.inc.static Plumed.cmake.static $(PLUMED_SCRIPTS_EXE)

lib-shared: $(PLUMED_KERNEL) $(PLUMED_SHARED_OBJ) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME) \
           Plumed.inc Plumed.inc.runtime Plumed.inc.shared \
           Plumed.cmake Plumed.cmake.runtime Plumed.cmake.shared

$(PLUMED_SCRIPTS_EXE):
	@echo "#! /bin/bash" > $@
	@echo "echo 'NOTE: shell only version, useful when plumed is cross compiled' >&2" >> $@
	@echo "export PLUMED_ROOT=\"\$${PLUMED_ROOT-$(PLUMED_ROOT)}\"" >> $@
	@echo "source \"\$$PLUMED_ROOT\"/scripts/$(subst plumed-,,$(@F)).sh" >> $@
	chmod a+x $@

$(PLUMED_KERNEL):$(OBJ_KERNEL)
# this workaround is needed to use a relative path in binding without
# overwriting the already installed library
# it moves the presently compiled version out of the way then move it back
# perhaps a cleaner solution is to use patchelf (linux) and install_name_tool (mac)
# but I am not sure how that would be portable
ifdef DESTDIR
	-test -f libplumedKernel.$(SOEXT) && mv libplumedKernel.$(SOEXT) .libplumedKernel.$(SOEXT)
	$(LDSO) -o libplumedKernel.$(SOEXT) $^ $(DYNAMIC_LIBS)
	mv libplumedKernel.$(SOEXT) $@
	-test -f .libplumedKernel.$(SOEXT) && mv .libplumedKernel.$(SOEXT) libplumedKernel.$(SOEXT)
else
	$(LDSO) -o $@ $^ $(DYNAMIC_LIBS)
endif

$(PLUMED_SHARED_OBJ): $(OBJ_KERNEL) $(OBJ_WRAPPER)
# this workaround is needed to use a relative path in binding without
# overwriting the already installed library
# it moves the presently compiled version out of the way then move it back
# perhaps a cleaner solution is to use patchelf (linux) and install_name_tool (mac)
# but I am not sure how that would be portable
ifdef DESTDIR
	-test -f libplumed.$(SOEXT) && mv libplumed.$(SOEXT) .libplumed.$(SOEXT)
	$(LDSO) -o libplumed.$(SOEXT) $^ $(DYNAMIC_LIBS)
	mv libplumed.$(SOEXT) $@
	-test -f .libplumed.$(SOEXT) && mv .libplumed.$(SOEXT) libplumed.$(SOEXT)
else
	$(LDSO) -o $@ $^ $(DYNAMIC_LIBS)
endif

$(PLUMED_MAIN_STATIC): $(OBJ_MAIN) $(OBJ_KERNEL) $(OBJ_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) $(DYNAMIC_LIBS)
	
$(PLUMED_MAIN_SHARED): $(OBJ_MAIN) $(PLUMED_SHARED_OBJ)
	$(LD) $(LDFLAGS) $(OBJ_MAIN) "$(realpath $(PLUMED_SHARED_OBJ))" -o $@ $(LIBS)

$(PLUMED_MAIN_RUNTIME): $(OBJ_MAIN) $(OBJ_DYNAMIC_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

Plumed.inc:
	@echo Building Plumed.inc
	@echo "PLUMED_RUNTIME_LOAD=" $(realpath $(OBJ_DYNAMIC_WRAPPER)) $(LIBS) $(LDFLAGS) > $@
	@echo "PLUMED_STATIC_LOAD=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER)) $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_SHARED_LOAD=" $(realpath .)/libplumed.$(SOEXT) $(LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_RUNTIME_DEPENDENCIES="   >> $@
	@echo "PLUMED_STATIC_DEPENDENCIES=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))  >> $@
	@echo "PLUMED_SHARED_DEPENDENCIES=" $(realpath .)/libplumed.$(SOEXT)  >> $@

Plumed.inc.runtime: Plumed.inc
	@echo Building Plumed.inc, runtime-linking version
	@echo "# PLUMED: runtime installation" > $@
	@grep PLUMED_RUNTIME $< | sed 's/PLUMED_RUNTIME_/PLUMED_/' >> $@

Plumed.inc.static: Plumed.inc
	@echo Building Plumed.inc, static version
	@echo "# PLUMED: static installation" > $@
	@grep PLUMED_STATIC $< | sed 's/PLUMED_STATIC_/PLUMED_/' >> $@

Plumed.inc.shared: Plumed.inc
	@echo Building Plumed.inc, shared version
	@echo "# PLUMED: shared installation" > $@
	@grep PLUMED_SHARED $< | sed 's/PLUMED_SHARED_/PLUMED_/' >> $@

Plumed.cmake:
	@echo Building Plumed.cmake
	@echo "set(PLUMED_RUNTIME_LOAD " $(realpath $(OBJ_DYNAMIC_WRAPPER)) $(LIBS) $(LDFLAGS)")" > $@
	@echo "set(PLUMED_STATIC_LOAD " $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER)) $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS)")" >> $@
	@echo "set(PLUMED_SHARED_LOAD " $(realpath ./)/libplumed.$(SOEXT) $(LIBS) $(LDFLAGS)")" >> $@
	@echo "set(PLUMED_RUNTIME_DEPENDENCIES) "   >> $@
	@echo "set(PLUMED_STATIC_DEPENDENCIES " $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))")"  >> $@
	@echo "set(PLUMED_SHARED_DEPENDENCIES " $(realpath ./)/libplumed.$(SOEXT)")"  >> $@

Plumed.cmake.runtime: Plumed.cmake
	@echo Building Plumed.cmake, runtime-linking version
	@echo "# PLUMED: runtime installation" > $@
	@grep PLUMED_RUNTIME $< | sed 's/PLUMED_RUNTIME_/PLUMED_/' >> $@

Plumed.cmake.static: Plumed.cmake
	@echo Building Plumed.cmake, static version
	@echo "# PLUMED: static installation" > $@
	@grep PLUMED_STATIC $< | sed 's/PLUMED_STATIC_/PLUMED_/' >> $@

Plumed.cmake.shared: Plumed.cmake
	@echo Building Plumed.cmake, shared version
	@echo "# PLUMED: shared installation" > $@
	@grep PLUMED_SHARED $< | sed 's/PLUMED_SHARED_/PLUMED_/' >> $@


# Individual dependencies

-include $(ALL_DEP)

# if machine dependent configuration has not been found:
else

.PHONY: error
error:
	@echo No configuration available
	@echo First run ./configure in the root directory
endif

# these targets are available anyway

.PHONY: clean dirslinks
clean:
	rm -f $(CLEANLIST)
	rm -f $(PLUMED_MAIN_STATIC) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME)
	rm -f Plumed.inc* plumed-* Plumed.cmake*
	rm -f .lib*

dirslinks:
	@echo "*** Updating links to used modules in all directories ***"
	$(MAKE) $(DIRSLINKS)
$(DIRSLINKS): %.links:
	$(MAKE) -C $* links


# generic makefile rules
include ../maketools/make.rules


