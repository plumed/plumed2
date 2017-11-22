
# version number
VERSION = $(shell if test -f ../../VERSION ; then grep -v \\\# ../../VERSION ; else echo "Unknown" ; fi )

# here we build the list of modules to be linked in plumedKernel
KERNEL_MODULES := $(addprefix ../,$(shell ../maketools/find_modules.sh))

# This adds modules which are built anyway and linked with some ad hoc procedure
DIRS := $(sort ../config ../wrapper ../main  $(KERNEL_MODULES))

# Targets to build links
DIRSLINKS:=$(addsuffix .links,$(DIRS))

# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../../Makefile.conf
endif

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

# when calling install-build, libraries in install should be built.
ifeq ($(PLUMED_INSTALL),Install)
  PLUMED_LIBDIR:=install/
else
  PLUMED_LIBDIR:=
endif

# plumed kernel
PLUMED_KERNEL := $(PLUMED_LIBDIR)libplumedKernel.$(SOEXT)

# plumed library (kernel + wrapper)
PLUMED_SHARED_OBJ := $(PLUMED_LIBDIR)libplumed.$(SOEXT)

# plumed wrapper library
ifdef AR_CR
PLUMED_WRAPPER_LIB := $(PLUMED_LIBDIR)libplumedWrapper.a
else
PLUMED_WRAPPER_LIB := 
endif

# plumed scripts
PLUMED_SCRIPTS_SH := $(shell cd ../../scripts/ ; ls *.sh 2>/dev/null)
PLUMED_SCRIPTS := $(PLUMED_SCRIPTS_SH:.sh=)
PLUMED_SCRIPTS_EXE := $(addprefix $(PLUMED_LIBDIR)plumed-,$(PLUMED_SCRIPTS))

# plumed main exec
PLUMED_MAIN_STATIC := $(PLUMED_LIBDIR)plumed
PLUMED_MAIN_SHARED := $(PLUMED_LIBDIR)plumed-shared
PLUMED_MAIN_RUNTIME := $(PLUMED_LIBDIR)plumed-runtime


.PHONY: all links $(DIRS) $(DIRSLINKS) all-static $(PLUMED_LIBDIR)Plumed.cmake $(PLUMED_LIBDIR)Plumed.inc install $(PLUMED_SCRIPTS_EXE) install-html
all:
	$(MAKE) dirslinks
	$(MAKE) dirs
ifdef SOEXT
	$(MAKE) lib
else
	$(MAKE) lib-static
endif
ifneq ($(PLUMED_INSTALL),Install)
	$(MAKE) install-build
endif


dirs:
	@echo "*** Compiling all directories ***"
	$(MAKE) $(DIRS)

install-build:
	mkdir -p install/
	$(MAKE) PLUMED_INSTALL=Install all
ifdef LD_RO
	$(MAKE) PLUMED_INSTALL=Install install/kernel.o
endif
# modulefile
	@cat modulefile.in | \
           sed "s|@_SOEXT_@|$(SOEXT)|" | \
           sed "s|@_program_name_@|$(program_name)|" | \
           sed "s|@_libdir_@|$(libdir)|" | \
           sed "s|@_bindir_@|$(bindir)|" | \
           sed "s|@_includedir_@|$(includedir)|" \
	   > install/modulefile
	@cat pkgconfig.in | \
           sed "s|@_SOEXT_@|$(SOEXT)|" | \
           sed "s|@_program_name_@|$(program_name)|" | \
           sed "s|@_libdir_@|$(libdir)|" | \
           sed "s|@_libs_@||" | \
           sed "s|@_bindir_@|$(bindir)|" | \
           sed "s|@_includedir_@|$(includedir)|" | \
           sed "s|@_prefix_@|$(prefix)|" | \
           sed "s|@_exec_prefix_@|$(exec_prefix)|" | \
           sed "s|@_VERSION_@|$(VERSION)|" | \
           sed "s|@_conflicts_@||" | \
           sed "s|@_cppflags_@||" | \
           sed "s|@_extra_@||" \
	   > install/pkgconfig
	@cat pkgconfig.in | \
           sed "s|@_SOEXT_@|$(SOEXT)|" | \
           sed "s|@_program_name_@|$(program_name)|" | \
           sed "s|@_libdir_@|$(libdir)|" | \
           sed "s|@_libs_@||" | \
           sed "s|@_bindir_@|$(bindir)|" | \
           sed "s|@_includedir_@|$(includedir)|" | \
           sed "s|@_prefix_@|$(prefix)|" | \
           sed "s|@_exec_prefix_@|$(exec_prefix)|" | \
           sed "s|@_VERSION_@|$(VERSION)|" | \
           sed "s|@_conflicts_@||" | \
           sed "s|@_cppflags_@|$(shell ../maketools/purge-d $(CPPFLAGS))|" | \
           sed "s|@_extra_@| (full version, including internals)|" \
	   > install/pkgconfigInternals
	@cat pkgconfig.in | \
           sed "s|@_SOEXT_@|$(SOEXT)|" | \
           sed "s|@_program_name_@|$(program_name)Wrapper|" | \
           sed "s|@_libdir_@|$(libdir)|" | \
           sed "s|@_libs_@|$(LIBS)|" | \
           sed "s|@_bindir_@|$(bindir)|" | \
           sed "s|@_includedir_@|$(includedir)|" | \
           sed "s|@_prefix_@|$(prefix)|" | \
           sed "s|@_exec_prefix_@|$(exec_prefix)|" | \
           sed "s|@_VERSION_@|$(VERSION)|" | \
           sed "s|@_conflicts_@||" | \
           sed "s|@_cppflags_@||" | \
           sed "s|@_extra_@| (wrapper-only version, set PLUMED_KERNEL at runtime) |" \
	   > install/pkgconfigWrapper

ifdef LD_RO
install/kernel.o: $(OBJ_KERNEL)
	$(LD_RO) install/kernel.o $(OBJ_KERNEL)
endif

# standard target (according to GNU doc)
install-html:
	if test -d ../../user-doc/html ; then mkdir -p "$(DESTDIR)$(htmldir)" && cd ../../ && tar cf - user-doc/html | tar xf - -C "$(DESTDIR)$(htmldir)/" ; fi
	if test -d ../../developer-doc/html ; then mkdir -p "$(DESTDIR)$(htmldir)" && cd ../../ && tar cf - developer-doc/html | tar xf - -C "$(DESTDIR)$(htmldir)/" ; fi
	if test -d ../../developer-doc/coverage ; then mkdir -p "$(DESTDIR)$(htmldir)" && cd ../../ && tar cf - developer-doc/coverage | tar xf - -C "$(DESTDIR)$(htmldir)/" ; fi
	if test -f ../../user-doc/manual.pdf ; then mkdir -p "$(DESTDIR)$(htmldir)/user-doc" && cp ../../user-doc/manual.pdf "$(DESTDIR)$(htmldir)/user-doc/" ; fi

install:
ifdef PLUMED_PREFIX
	@echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	@echo "WARNING: using PLUMED_PREFIX variable is not allowed anymore"
	@echo "         please use one of the following choices:"
	@echo "         (1) at configure time:"
	@echo "           ./configure --prefix=$(PLUMED_PREFIX)"
	@echo "         (2) or later, at install time:"
	@echo "           make install prefix=$(PLUMED_PREFIX)"
	@echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
endif
ifdef PLUMED_LIBSUFFIX
	@echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	@echo "WARNING: using PLUMED_LIBSUFFIX variable is not allowed anymore"
	@echo "         please use the following option at configure time"
	@echo "         ./configure --program-suffix='_$(PLUMED_LIBSUFFIX)'"
	@echo "         there is no way to change the suffix at install time"
	@echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
endif
ifndef program_name
	$(error "it seems Makefile.conf does not defined program_name - please configure again!")
endif
ifndef bindir
	$(error "it seems Makefile.conf does not defined bindir - please configure again!")
endif
ifndef libdir
	$(error "it seems Makefile.conf does not defined libdir - please configure again!")
endif
ifndef includedir
	$(error "it seems Makefile.conf does not defined includedir - please configure again!")
endif
ifndef htmldir
	$(error "it seems Makefile.conf does not defined htmldir - please configure again!")
endif
	cd ../config ; $(MAKE)
# in principle we should uninstall before install
# however, I think it is safer not to do it here since users could make mistake
# and remove important files by just typing "sudo make install" with a 
# malformed prefix variables
	# $(MAKE) uninstall
	mkdir -p "$(DESTDIR)$(bindir)"
	mkdir -p "$(DESTDIR)$(libdir)"
	mkdir -p "$(DESTDIR)$(libdir)/pkgconfig"
	mkdir -p "$(DESTDIR)$(includedir)/$(program_name)"
	mkdir -p "$(DESTDIR)$(libdir)/$(program_name)/obj"
	mkdir -p "$(DESTDIR)$(libdir)/$(program_name)/src/lib"
	mkdir -p "$(DESTDIR)$(libdir)/$(program_name)/src/config"
	mkdir -p "$(DESTDIR)$(libdir)/$(program_name)/src/colvar"
# rebuild libraries and main using new target dir
	$(MAKE) install-build
# copy documentation:
	$(MAKE) install-html
# copy patches and scripts
	cd ../../ ; tar cf - patches scripts | tar xf - -C "$(DESTDIR)$(libdir)/$(program_name)/"
# copy Template.cpp file (we leave it in src/colvar/ for backward compatibility)
	cp ../colvar/Template.cpp "$(DESTDIR)$(libdir)/$(program_name)/src/colvar/"
# copy compile_options.sh file (we leave it in src/colvar/ for backward compatibility)
	cp ../config/compile_options.sh "$(DESTDIR)$(libdir)/$(program_name)/src/config/compile_options.sh"
# copy config.txt file (we leave it in src/colvar/ for backward compatibility)
	cp ../config/config.txt "$(DESTDIR)$(libdir)/$(program_name)/src/config/config.txt"
ifdef LD_RO
	cp install/kernel.o "$(DESTDIR)$(libdir)/$(program_name)/obj/kernel.o"
else
	../maketools/copyobjects "$(DESTDIR)$(libdir)/$(program_name)/obj/k" $(subst ../config/Config.o,../config/ConfigInstall.o,$(OBJ_KERNEL))
endif
	cp $(OBJ_WRAPPER) "$(DESTDIR)$(libdir)/$(program_name)/obj/PlumedStatic.o"
	cp $(OBJ_DYNAMIC_WRAPPER) "$(DESTDIR)$(libdir)/$(program_name)/obj/Plumed.o"
# also copy .h files into include/ dir
	cd ../../src ; tar cf - */*.h */*/*.h | tar xf - -C "$(DESTDIR)$(includedir)/$(program_name)/"
# install executable
	cp install/plumed "$(DESTDIR)$(bindir)/$(program_name)"
# patch file for cross compiling
	cp install/plumed-patch "$(DESTDIR)$(bindir)/$(program_name)-patch"
# config file for cross compiling
	cp install/plumed-config "$(DESTDIR)$(bindir)/$(program_name)-config"
# copy include files for cmake and make
# we leave it in src/lib/ for backward compatibility
	cp install/Plumed.cmake* "$(DESTDIR)$(libdir)/$(program_name)/src/lib/"
	cp install/Plumed.inc* "$(DESTDIR)$(libdir)/$(program_name)/src/lib/"
# other shell-based cltools are copied
	for file in $(PLUMED_SCRIPTS) ; do cp install/plumed-$$file "$(DESTDIR)$(libdir)/$(program_name)/plumed-$$file" ; done
ifdef SOEXT
	cp install/libplumed.$(SOEXT) "$(DESTDIR)$(libdir)/lib$(program_name).$(SOEXT)"
	cp install/libplumedKernel.$(SOEXT) "$(DESTDIR)$(libdir)/lib$(program_name)Kernel.$(SOEXT)"
ifeq ($(SOEXT),dylib)
	install_name_tool -id "$(libdir)/lib$(program_name).$(SOEXT)" "$(DESTDIR)$(libdir)/lib$(program_name).$(SOEXT)"
	install_name_tool -id "$(libdir)/lib$(program_name)Kernel.$(SOEXT)" "$(DESTDIR)$(libdir)/lib$(program_name)Kernel.$(SOEXT)"
endif
endif
ifdef AR_CR
	cp install/libplumedWrapper.a "$(DESTDIR)$(libdir)/lib$(program_name)Wrapper.a"
endif
# modulefile
	cp install/modulefile "$(DESTDIR)$(libdir)/$(program_name)/modulefile"
# pkgconfig files
	cp install/pkgconfig "$(DESTDIR)$(libdir)/pkgconfig/$(program_name).pc"
ifdef AR_CR
	cp install/pkgconfigWrapper "$(DESTDIR)$(libdir)/pkgconfig/$(program_name)Wrapper.pc"
endif
	cp install/pkgconfigInternals "$(DESTDIR)$(libdir)/pkgconfig/$(program_name)Internals.pc"
# vim syntax file
	if test -d ../../vim/syntax ; then mkdir -p "$(DESTDIR)$(libdir)/$(program_name)" && cd ../../ && tar cf - vim/syntax | tar xf - -C "$(DESTDIR)$(libdir)/$(program_name)/" ; fi
	if test -d ../../vim/help ; then mkdir -p "$(DESTDIR)$(libdir)/$(program_name)" && cd ../../ && tar cf - vim/help | tar xf - -C "$(DESTDIR)$(libdir)/$(program_name)/" ; fi
	if test -f ../../vim/scripts.vim ; then cp ../../vim/scripts.vim  "$(DESTDIR)$(libdir)/$(program_name)/vim/" ; fi
# making everything visible:
	chmod -R go+rX,go-w "$(DESTDIR)$(libdir)/$(program_name)"
	chmod -R go+rX,go-w "$(DESTDIR)$(includedir)/$(program_name)"
	chmod -R go+rX,go-w "$(DESTDIR)$(bindir)/$(program_name)"
	chmod -R go+rX,go-w "$(DESTDIR)$(bindir)/$(program_name)-patch"
	@echo
	@echo "*** PLUMED has been installed ***"
	@echo
ifdef DESTDIR
	@echo "WARNING: Staged install with DESTDIR = $(DESTDIR)"
endif
	@echo "Install prefix : $(prefix)"
	@echo "Full name      : $(program_name)"
	@echo
	@echo "Setup your environment"
	@echo "- Ensure this is in your execution path         : $(bindir)"
	@echo "- Ensure this is in your include path           : $(includedir)"
ifdef SOEXT
	@echo "- Ensure this is in your library path           : $(libdir)"
	@echo "- Ensure this is in your PKG_CONFIG_PATH path   : $(libdir)/pkgconfig"
	@echo "For runtime binding:"
	@echo "- Set this environment variable                 : PLUMED_KERNEL=$(libdir)/lib$(program_name)Kernel.$(SOEXT)"
endif
	@echo
	@echo "To create a tcl module that sets all the variables above, use this one as a starting point:"
	@echo "$(libdir)/$(program_name)/modulefile"
	@echo
	@echo "To uninstall, remove the following files and directories:"
	@echo $(libdir)/$(program_name)
	@echo $(htmldir)
	@echo $(includedir)/$(program_name)
	@echo $(bindir)/$(program_name)
	@echo $(bindir)/$(program_name)-patch
	@echo $(bindir)/$(program_name)-config
	@echo $(libdir)/pkgconfig/$(program_name).pc
ifdef SOEXT
	@echo $(libdir)/lib$(program_name).$(SOEXT)
	@echo $(libdir)/lib$(program_name)Kernel.$(SOEXT)
endif
	@if test -d ../../vim/help ; then echo "A vim plugin can be found here: $(libdir)/$(program_name)/vim/" ; fi
	@if test -d ../../vim/help ; then echo "Copy it to $(HOME)/.vim/ directory" ; fi
	@if test -d ../../vim/help ; then echo "Alternatively:" ; fi
	@if test -d ../../vim/help ; then echo "- Set this environment variable         : PLUMED_VIMPATH=$(libdir)/$(program_name)/vim" ; fi
	@if test -d ../../vim/help ; then echo "- Add the command 'let &runtimepath.=','.\$$PLUMED_VIMPATH' to your .vimrc file" ; fi
	@if test -d ../../vim/help ; then echo "From vim, you can use :set syntax=$(program_name) to enable it" ; fi
ifeq ($(program_can_run),no)
	@echo "WARNING: $(program_name) executable will not run on this machine"
	@echo "WARNING: to patch an MD code use 'plumed-patch'"
	@echo "WARNING: to verify configuration use 'plumed-config'"
	@echo "WARNING: a few other command line tools are avalable as $(libdir)/$(program_name)/plumed-*"
else ifeq ($(program_can_run_mpi),no)
	@echo "WARNING: $(program_name) executable will not run on this machine"
	@echo "WARNING: unless you invoke it as '$(program_name) --no-mpi'"
	@echo "WARNING: to patch an MD code now use '$(program_name) --no-mpi patch'"
	@echo "WARNING: (notice that MPI will be available anyway in the patched code)"
	@echo "WARNING: all command line tools are available as '$(program_name) --no-mpi name-of-the-tool'"
	@echo "WARNING: e.g. '$(program_name) --no-mpi driver'"
endif


uninstall:
	@echo WARNING
	@echo I will remove the following files now
	@echo "$(DESTDIR)$(libdir)/$(program_name)"
	@echo "$(DESTDIR)$(htmldir)"
	@echo "$(DESTDIR)$(includedir)/$(program_name)"
	@echo "$(DESTDIR)$(bindir)/$(program_name)"
	@echo "$(DESTDIR)$(bindir)/$(program_name)-patch"
	@echo "$(DESTDIR)$(bindir)/$(program_name)-config"
	@echo "$(DESTDIR)$(libdir)/lib$(program_name).$(SOEXT)"
	@echo "$(DESTDIR)$(libdir)/lib$(program_name)Kernel.$(SOEXT)"
	../maketools/confirm_uninstall
	@echo Removing previously installed files
	rm -fr "$(DESTDIR)$(libdir)/$(program_name)"
	rm -fr "$(DESTDIR)$(htmldir)"
	rm -fr "$(DESTDIR)$(includedir)/$(program_name)"
	rm -f "$(DESTDIR)$(bindir)/$(program_name)"
	rm -f "$(DESTDIR)$(bindir)/$(program_name)-patch"
	rm -f "$(DESTDIR)$(bindir)/bin/$(program_name)-config"
	rm -f "$(DESTDIR)$(libdir)/lib$(program_name).$(SOEXT)"
	rm -f "$(DESTDIR)$(libdir)/lib$(program_name)Kernel.$(SOEXT)"
	rm -f "$(DESTDIR)$(libdir)/lib/pkgconfig/$(PLUMED_PROGRAM_NAME).pc"

$(DIRS):
	$(MAKE) -C $@ obj

.PHONY: lib lib-static lib-shared lib-install lib-static-install

lib: lib-static lib-shared

lib-static: $(PLUMED_MAIN_STATIC) $(OBJ_WRAPPER) $(PLUMED_LIBDIR)Plumed.inc.static $(PLUMED_LIBDIR)Plumed.cmake.static $(PLUMED_SCRIPTS_EXE)

lib-shared: $(PLUMED_KERNEL) $(PLUMED_SHARED_OBJ) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME) $(PLUMED_WRAPPER_LIB) \
           $(PLUMED_LIBDIR)Plumed.inc $(PLUMED_LIBDIR)Plumed.inc.runtime $(PLUMED_LIBDIR)Plumed.inc.shared \
           $(PLUMED_LIBDIR)Plumed.cmake $(PLUMED_LIBDIR)Plumed.cmake.runtime $(PLUMED_LIBDIR)Plumed.cmake.shared

$(PLUMED_SCRIPTS_EXE):
	@echo "#! /bin/bash" > $@
	@echo "echo 'NOTE: shell only version, useful when plumed is cross compiled' >&2" >> $@
# The logic implemented here should be the same as in config/Config.inc.in
# Please keep them consistent
ifeq ($(PLUMED_INSTALL),Install)
# If plumed is installed, single paths can be overridden using env vars
	@echo "export PLUMED_IS_INSTALLED=yes" >> $@
	@echo "export PLUMED_ROOT=\"\$${PLUMED_ROOT-$(libdir)/$(program_name)}\"" >> $@
	@echo "export PLUMED_INCLUDEDIR=\"\$${PLUMED_INCLUDEDIR-$(includedir)}\"" >> $@
	@echo "export PLUMED_HTMLDIR=\"\$${PLUMED_HTMLDIR=$(htmldir)}\"" >> $@
	@echo "export PLUMED_PROGRAM_NAME=\"\$${PLUMED_PROGRAM_NAME-$(program_name)}\"" >> $@
	@echo "source \"\$$PLUMED_ROOT\"/scripts/$(subst plumed-,,$(@F)).sh" >> $@
else
	@echo "export PLUMED_IS_INSTALLED=no" >> $@
	@echo "export PLUMED_ROOT=\"\$${PLUMED_ROOT-$(realpath ../../)}\"" >> $@
	@echo "export PLUMED_INCLUDEDIR=\"\$${PLUMED_ROOT}/include\"" >> $@
	@echo "export PLUMED_HTMLDIR=\"\$${PLUMED_ROOT}\"" >> $@
	@echo "export PLUMED_PROGRAM_NAME=plumed" >> $@
	@echo "source \"\$$PLUMED_ROOT\"/scripts/$(subst plumed-,,$(@F)).sh" >> $@
endif
	chmod a+x $@

$(PLUMED_KERNEL):$(OBJ_KERNEL)
ifeq ($(PLUMED_INSTALL),Install)
	mkdir -p install
	cd install && $(LDSHARED) -o libplumedKernel.$(SOEXT) $(addprefix ../,$(OBJ_KERNEL)) $(DYNAMIC_LIBS)
else
	$(LDSHARED) -o $@ $^ $(DYNAMIC_LIBS)
ifeq ($(SOEXT),dylib)
	install_name_tool -id "$(realpath .)/libplumedKernel.$(SOEXT)" libplumedKernel.$(SOEXT)
endif
endif

$(PLUMED_SHARED_OBJ): $(OBJ_KERNEL) $(OBJ_WRAPPER)
ifeq ($(PLUMED_INSTALL),Install)
	mkdir -p install
	cd install && $(LDSHARED) -o libplumed.$(SOEXT) $(addprefix ../,$(OBJ_KERNEL)) $(addprefix ../,$(OBJ_WRAPPER)) $(DYNAMIC_LIBS)
else
	$(LDSHARED) -o $@ $^ $(DYNAMIC_LIBS)
ifeq ($(SOEXT),dylib)
	install_name_tool -id "$(realpath .)/libplumed.$(SOEXT)" $@
endif
endif

ifdef AR_CR
$(PLUMED_WRAPPER_LIB): $(OBJ_DYNAMIC_WRAPPER)
ifeq ($(PLUMED_INSTALL),Install)
	mkdir -p install
	cd install && $(AR_CR) libplumedWrapper.a $(addprefix ../,$(OBJ_DYNAMIC_WRAPPER))
else
	$(AR_CR) $@ $(OBJ_DYNAMIC_WRAPPER)
endif
endif

$(PLUMED_MAIN_STATIC): $(OBJ_MAIN) $(OBJ_KERNEL) $(OBJ_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) $(DYNAMIC_LIBS)
	
$(PLUMED_MAIN_SHARED): $(OBJ_MAIN) $(PLUMED_SHARED_OBJ)
	$(LD) $(LDFLAGS) $(OBJ_MAIN) "$(realpath $(PLUMED_SHARED_OBJ))" -o $@ $(LIBS)

$(PLUMED_MAIN_RUNTIME): $(OBJ_MAIN) $(OBJ_DYNAMIC_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

ifeq ($(PLUMED_INSTALL),Install)

$(PLUMED_LIBDIR)Plumed.inc:
	@echo Building Plumed.inc
ifdef AR_CR
	@echo "PLUMED_RUNTIME_LOAD= \"$(libdir)/lib$(program_name)Wrapper.a\" $(LIBS) $(LDFLAGS)" > $@
else
	@echo "PLUMED_RUNTIME_LOAD= \"$(libdir)/$(program_name)/obj/Plumed.o\" $(LIBS) $(LDFLAGS)" > $@
endif
ifdef LD_RO
	@echo "PLUMED_STATIC_LOAD= \"$(libdir)/$(program_name)/obj/kernel.o\" \"$(libdir)/$(program_name)/obj/PlumedStatic.o\" $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS)" >> $@
else
# single quote required to preserve double quote in resulting file
	@echo "PLUMED_STATIC_LOAD= '$(shell ../maketools/listobjects "$(libdir)/$(program_name)/obj/k" $(OBJ_KERNEL))' \"$(libdir)/$(program_name)/obj/PlumedStatic.o\" $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS)" >> $@
endif
	@echo "PLUMED_SHARED_LOAD= \"$(libdir)/lib$(program_name).$(SOEXT)\" $(LIBS) $(LDFLAGS)" >> $@
	@echo "PLUMED_RUNTIME_DEPENDENCIES="   >> $@
# in principle all objects are replaced at the same time
# however, to keep compatibility with the include files generated in plumed 2.0 and 2.1,
# I put all the objects (or shared objects) here
ifdef LD_RO
	@echo "PLUMED_STATIC_DEPENDENCIES= \"$(libdir)/$(program_name)/obj/kernel.o\" \"$(libdir)/$(program_name)/obj/PlumedStatic.o\"" >> $@
else
# single quote required to preserve double quote in resulting file
	@echo "PLUMED_STATIC_DEPENDENCIES= '$(shell ../maketools/listobjects "$(libdir)/$(program_name)/obj/k" $(OBJ_KERNEL))' \"$(libdir)/$(program_name)/obj/PlumedStatic.o\"" >> $@
endif
	@echo "PLUMED_SHARED_DEPENDENCIES= \"$(libdir)/lib$(program_name).$(SOEXT)\""  >> $@
else
Plumed.inc:
	@echo Building Plumed.inc
ifdef AR_CR
	@echo "PLUMED_RUNTIME_LOAD=" $(realpath .)/libplumedWrapper.a $(LIBS) $(LDFLAGS) > $@
else
	@echo "PLUMED_RUNTIME_LOAD=" $(realpath $(OBJ_DYNAMIC_WRAPPER)) $(LIBS) $(LDFLAGS) > $@
endif
	@echo "PLUMED_STATIC_LOAD=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER)) $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_SHARED_LOAD=" $(realpath .)/libplumed.$(SOEXT) $(LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_RUNTIME_DEPENDENCIES="   >> $@
	@echo "PLUMED_STATIC_DEPENDENCIES=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))  >> $@
	@echo "PLUMED_SHARED_DEPENDENCIES=" $(realpath .)/libplumed.$(SOEXT)  >> $@
endif

$(PLUMED_LIBDIR)Plumed.cmake: $(PLUMED_LIBDIR)Plumed.inc
	@echo Building Plumed.cmake
	@../maketools/inc2cmake < $(PLUMED_LIBDIR)Plumed.inc > $@

$(PLUMED_LIBDIR)Plumed.inc.runtime: $(PLUMED_LIBDIR)Plumed.inc
	@echo Building Plumed.inc, runtime-linking version
	@echo "# PLUMED: runtime installation" > $@
	@grep PLUMED_RUNTIME $< | sed 's/PLUMED_RUNTIME_/PLUMED_/' >> $@

$(PLUMED_LIBDIR)Plumed.inc.static: $(PLUMED_LIBDIR)Plumed.inc
	@echo Building Plumed.inc, static version
	@echo "# PLUMED: static installation" > $@
	@grep PLUMED_STATIC $< | sed 's/PLUMED_STATIC_/PLUMED_/' >> $@

$(PLUMED_LIBDIR)Plumed.inc.shared: $(PLUMED_LIBDIR)Plumed.inc
	@echo Building Plumed.inc, shared version
	@echo "# PLUMED: shared installation" > $@
	@grep PLUMED_SHARED $< | sed 's/PLUMED_SHARED_/PLUMED_/' >> $@

$(PLUMED_LIBDIR)Plumed.cmake.runtime: $(PLUMED_LIBDIR)Plumed.cmake
	@echo Building Plumed.cmake, runtime-linking version
	@echo "# PLUMED: runtime installation" > $@
	@grep PLUMED_RUNTIME $< | sed 's/PLUMED_RUNTIME_/PLUMED_/' >> $@

$(PLUMED_LIBDIR)Plumed.cmake.static: $(PLUMED_LIBDIR)Plumed.cmake
	@echo Building Plumed.cmake, static version
	@echo "# PLUMED: static installation" > $@
	@grep PLUMED_STATIC $< | sed 's/PLUMED_STATIC_/PLUMED_/' >> $@

$(PLUMED_LIBDIR)Plumed.cmake.shared: $(PLUMED_LIBDIR)Plumed.cmake
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
# this is not necessary, I leave it here to remove spurious plumed.vim files generated with
# the first version of plumed.vim syntax file:
	rm -f plumed.vim
	rm -f .lib*
	rm -fr install/

dirslinks:
	@echo "*** Updating links to used modules in all directories ***"
	$(MAKE) $(DIRSLINKS)
$(DIRSLINKS): %.links:
	$(MAKE) -C $* links


# generic makefile rules
include ../maketools/make.rules


