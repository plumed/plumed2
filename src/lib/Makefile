# here we build the list of modules to be linked in plumedKernel
KERNEL_MODULES=$(addprefix ../,$(shell ../maketools/find_modules.sh))

# This adds modules which are built anyway and linked with some ad hoc procedure
DIRS=$(sort ../config ../wrapper ../main  $(KERNEL_MODULES))


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
SRC_KERNEL=$(shell for dir in $(KERNEL_MODULES) ; do ls $$dir/*.cpp 2>/dev/null ; ls $$dir/*.c 2>/dev/null ; done) ../config/Config$(PLUMED_INSTALL).o

# corresponding objects:
OBJ_WRAPPER=../wrapper/PlumedStatic.o
OBJ_DYNAMIC_WRAPPER=../wrapper/Plumed.o
OBJ_MAIN=$(SRC_MAIN:.cpp=.o)
XXX=$(SRC_KERNEL:.cpp=.o)
OBJ_KERNEL=$(XXX:.c=.o)
ALL_OBJ=$(OBJ_KERNEL) $(OBJ_WRAPPER) $(OBJ_DYNAMIC_WRAPPER) $(OBJ_MAIN)

# plumed kernel
PLUMED_KERNEL=$(PLUMED_LIBDIR)libplumedKernel.$(SOEXT)

# plumed library (kernel + wrapper)
PLUMED_SHARED_OBJ=$(PLUMED_LIBDIR)libplumed.$(SOEXT)

# plumed main exec
PLUMED_MAIN_STATIC=$(PLUMED_LIBDIR)plumed
PLUMED_MAIN_SHARED=$(PLUMED_LIBDIR)plumed-shared
PLUMED_MAIN_RUNTIME=$(PLUMED_LIBDIR)plumed-runtime
PLUMED_MAIN_PATCH=$(PLUMED_LIBDIR)plumed-patch

.PHONY: all dirs $(DIRS) all-static
ifdef SOEXT
all:
	@echo First making other dirs $(DIRS)
	$(MAKE) dirs
	$(MAKE) lib
else
all:
	$(MAKE) all-static
endif

all-static:
	@echo First making other dirs $(DIRS)
	$(MAKE) dirs
	$(MAKE) lib-static


install:
	cd ../config ; $(MAKE)
# two steps to be sure that ../config/install.conf is up to date
	$(MAKE) install-do

install-do:
	$(MAKE) uninstall
	mkdir -p $(PLUMED_INSTALL_PREFIX)/bin
	mkdir -p $(PLUMED_INSTALL_ROOT)/src/lib
	mkdir -p $(PLUMED_INSTALL_PREFIX)/include/plumed$(PLUMED_INSTALL_LIBSUFFIX)
# rebuild libraries and main using new target dir
	$(MAKE) PLUMED_LIBDIR=$(PLUMED_INSTALL_ROOT)/src/lib/ PLUMED_INSTALL=Install all
# copy .h, .cpp .o compile_options.sh and documentation:
	cd ../../ ; tar cf - patches scripts src/*/*.h src/*/*.cpp src/*/*.o src/config/compile_options.sh *-doc | tar xf - -C $(PLUMED_INSTALL_ROOT)/
# also copy .h files into include/ dir
	cd ../../src ; tar cf - */*.h $$(find . -type l) | tar xf - -C $(PLUMED_INSTALL_PREFIX)/include/plumed$(PLUMED_INSTALL_LIBSUFFIX)/
# copy inc files changing the path
	for file in Plumed.inc* Plumed.cmake* ; do sed "s| $(realpath ../../)| $(PLUMED_INSTALL_ROOT)/|g" $$file > $(PLUMED_INSTALL_ROOT)/src/lib/$$file ; done
# modulefile
	@echo "#%Module1.0##############################################"          >  $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "proc ModulesHelp { } {"                                             >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "  puts stderr \"\tAdds paths for Plumed\""                          >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "}"                                                                  >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "module-whatis   \"Adds paths for Plumed\""                          >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "# Manually add here dependencies and conflicts"                     >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "prepend-path  PATH               $(PLUMED_INSTALL_PREFIX)/bin"      >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "prepend-path  CPATH              $(PLUMED_INSTALL_PREFIX)/include/" >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "prepend-path  INCLUDE            $(PLUMED_INSTALL_PREFIX)/include/" >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "prepend-path  LIBRARY_PATH       $(PLUMED_INSTALL_PREFIX)/lib/"     >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "prepend-path  LD_LIBRARY_PATH    $(PLUMED_INSTALL_PREFIX)/lib/"     >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
	@echo "prepend-path  DYLD_LIBRARY_PATH  $(PLUMED_INSTALL_PREFIX)/lib/"     >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
ifdef SOEXT
	@echo "setenv        PLUMED_KERNEL      $(PLUMED_INSTALL_PREFIX)/lib/libplumedKernel.$(SOEXT)" >> $(PLUMED_INSTALL_ROOT)/src/lib/modulefile
endif
# install links
	ln -s $(PLUMED_INSTALL_ROOT)/src/lib/plumed $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)
	ln -s $(PLUMED_INSTALL_ROOT)/src/lib/plumed-patch $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)-patch
ifdef SOEXT
	ln -s $(PLUMED_INSTALL_ROOT)/src/lib/libplumed.$(SOEXT) $(PLUMED_INSTALL_PREFIX)/lib/libplumed$(PLUMED_INSTALL_LIBSUFFIX).$(SOEXT)
	ln -s $(PLUMED_INSTALL_ROOT)/src/lib/libplumedKernel.$(SOEXT) $(PLUMED_INSTALL_PREFIX)/lib/libplumed$(PLUMED_INSTALL_LIBSUFFIX)Kernel.$(SOEXT)
endif
	chmod -R go+rX,go-w $(PLUMED_INSTALL_ROOT)
	chmod -R go+rX,go-w $(PLUMED_INSTALL_PREFIX)/include/plumed$(PLUMED_INSTALL_LIBSUFFIX)
	chmod -R go+rX,go-w $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)
	chmod -R go+rX,go-w $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)-patch
	@echo "To uninstall, remove the following files and directories:"
	@echo  $(PLUMED_INSTALL_ROOT)
	@echo  $(PLUMED_INSTALL_PREFIX)/include/plumed$(PLUMED_INSTALL_LIBSUFFIX)
	@echo $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)
	@echo $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)-patch
ifdef SOEXT
	@echo $(PLUMED_INSTALL_PREFIX)/lib/libplumed$(PLUMED_INSTALL_LIBSUFFIX).$(SOEXT)
	@echo $(PLUMED_INSTALL_PREFIX)/lib/libplumed$(PLUMED_INSTALL_LIBSUFFIX)Kernel.$(SOEXT)
endif

uninstall:
	rm -fr $(PLUMED_INSTALL_ROOT)
	rm -fr $(PLUMED_INSTALL_PREFIX)/include/plumed$(PLUMED_INSTALL_LIBSUFFIX)
	rm -f $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)
	rm -f $(PLUMED_INSTALL_PREFIX)/bin/plumed$(PLUMED_INSTALL_LIBSUFFIX)-patch
ifdef SOEXT
	rm -f $(PLUMED_INSTALL_PREFIX)/lib/libplumed$(PLUMED_INSTALL_LIBSUFFIX).$(SOEXT)
	rm -f $(PLUMED_INSTALL_PREFIX)/lib/libplumed$(PLUMED_INSTALL_LIBSUFFIX)Kernel.$(SOEXT)
endif

dirs: $(DIRS)

$(DIRS):
	$(MAKE) -C $@


# default target
.PHONY: lib lib-static lib-shared

lib: lib-static lib-shared

lib-static: $(PLUMED_MAIN_STATIC) $(OBJ_WRAPPER) Plumed.inc.static Plumed.cmake.static $(PLUMED_MAIN_PATCH)

lib-shared: $(PLUMED_KERNEL) $(PLUMED_SHARED_OBJ) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME) \
           Plumed.inc Plumed.inc.runtime Plumed.inc.shared \
           Plumed.cmake Plumed.cmake.runtime Plumed.cmake.shared

$(PLUMED_MAIN_PATCH):
	@echo "#! /bin/bash" > $(PLUMED_MAIN_PATCH)
	@echo "echo 'NOTE: shell only version of path, useful for cross compiling'" >> $(PLUMED_MAIN_PATCH)
	@echo "export PLUMED_ROOT=\"$(realpath $(PLUMED_LIBDIR)../../)\"" >> $(PLUMED_MAIN_PATCH)
	@echo "source \"\$$PLUMED_ROOT\"/scripts/patch.sh" >> $(PLUMED_MAIN_PATCH)
	chmod a+x $(PLUMED_MAIN_PATCH)

$(PLUMED_KERNEL):$(OBJ_KERNEL)
	$(LDSO) -o $@ $^ $(DYNAMIC_LIBS)

$(PLUMED_SHARED_OBJ): $(OBJ_KERNEL) $(OBJ_WRAPPER)
	$(LDSO) -o $@ $^ $(DYNAMIC_LIBS)

$(PLUMED_MAIN_STATIC): $(OBJ_MAIN) $(OBJ_KERNEL) $(OBJ_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) $(DYNAMIC_LIBS)
	
$(PLUMED_MAIN_SHARED): $(OBJ_MAIN) $(PLUMED_SHARED_OBJ)
	$(LD) $(LDFLAGS) $(OBJ_MAIN) $(realpath $(PLUMED_SHARED_OBJ)) -o $@ $(LIBS)

$(PLUMED_MAIN_RUNTIME): $(OBJ_MAIN) $(OBJ_DYNAMIC_WRAPPER)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

Plumed.inc:
	@echo Building Plumed.inc
	@echo "PLUMED_RUNTIME_LOAD=" $(realpath $(OBJ_DYNAMIC_WRAPPER)) $(LIBS) $(LDFLAGS) > $@
	@echo "PLUMED_STATIC_LOAD=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER)) $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_SHARED_LOAD=" $(realpath ./)/$(PLUMED_SHARED_OBJ) $(LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_RUNTIME_DEPENDENCIES="   >> $@
	@echo "PLUMED_STATIC_DEPENDENCIES=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))  >> $@
	@echo "PLUMED_SHARED_DEPENDENCIES=" $(realpath ./)/$(PLUMED_SHARED_OBJ)  >> $@

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
	@echo "set(PLUMED_SHARED_LOAD " $(realpath ./)/$(PLUMED_SHARED_OBJ) $(LIBS) $(LDFLAGS)")" >> $@
	@echo "set(PLUMED_RUNTIME_DEPENDENCIES) "   >> $@
	@echo "set(PLUMED_STATIC_DEPENDENCIES " $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))")"  >> $@
	@echo "set(PLUMED_SHARED_DEPENDENCIES " $(realpath ./)/$(PLUMED_SHARED_OBJ)")"  >> $@

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
	@echo First run ./configure.sh in the root directory
endif

# this target is available anyway

.PHONY: clean
clean:
	rm -f $(CLEANLIST)
	rm -f $(PLUMED_MAIN_STATIC) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME)
	rm -f Plumed.inc* plumed-* Plumed.cmake*

# generic makefile rules
include ../maketools/make.rules


