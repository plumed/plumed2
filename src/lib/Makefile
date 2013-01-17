# here we build the list of modules to be linked in plumedKernel
KERNEL_MODULES=$(addprefix ../,$(shell ../maketools/find_modules.sh))

# This adds modules which are built anyway and linked with some ad hoc procedure
DIRS=$(sort ../config ../wrapper ../main  $(KERNEL_MODULES))


# include the machine dependent configuration
-include ../../Makefile.conf

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
SRC_KERNEL=$(shell for dir in $(KERNEL_MODULES) ; do ls $$dir/*.{c,cpp} 2>/dev/null ; done) ../config/Config.o

# corresponding objects:
OBJ_WRAPPER=../wrapper/PlumedStatic.o
OBJ_DYNAMIC_WRAPPER=../wrapper/Plumed.o
OBJ_MAIN=$(SRC_MAIN:.cpp=.o)
XXX=$(SRC_KERNEL:.cpp=.o)
OBJ_KERNEL=$(XXX:.c=.o)
ALL_OBJ=$(OBJ_KERNEL) $(OBJ_WRAPPER) $(OBJ_DYNAMIC_WRAPPER) $(OBJ_MAIN)

# plumed kernel
PLUMED_KERNEL=libplumedKernel.$(SOEXT)

# plumed library (kernel + wrapper)
PLUMED_SHARED_OBJ=libplumed.$(SOEXT)

# plumed main exec
PLUMED_MAIN_STATIC=plumed
PLUMED_MAIN_SHARED=plumed-shared
PLUMED_MAIN_RUNTIME=plumed-runtime

.PHONY: all dirs $(DIRS)
all:
	@echo First making other dirs $(DIRS)
	$(MAKE) dirs
	$(MAKE) lib

dirs: $(DIRS)

$(DIRS):
	$(MAKE) -C $@


# default target
.PHONY: lib
lib: $(PLUMED_KERNEL) $(PLUMED_SHARED_OBJ) $(OBJ_WRAPPER) \
     $(PLUMED_MAIN_STATIC) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME) \
     Plumed.inc Plumed.inc.runtime Plumed.inc.shared Plumed.inc.static \
     Plumed.cmake Plumed.cmake.runtime Plumed.cmake.shared Plumed.cmake.static \
     plumed-patch

plumed-patch:
	@echo "#! /bin/bash" > plumed-patch
	@echo "echo 'NOTE: shell only version of path, useful for cross compiling'" >> plumed-patch
	@echo "export PLUMED_ROOT=\"$(realpath ../../)\"" >> plumed-patch
	@echo "source \"$(realpath ../../scripts/patch.sh)\"" >> plumed-patch
	chmod a+x plumed-patch

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

Plumed.inc: $(OBJ_KERNEL) $(OBJ_WRAPPER) $(PLUMED_SHARED_OBJ) $(OBJ_DYNAMIC_WRAPPER)
	@echo Building Plumed.inc
	@echo "PLUMED_INCLUDE="-I$(realpath ../../include) > $@
	@echo "PLUMED_RUNTIME_LOAD=" $(realpath $(OBJ_DYNAMIC_WRAPPER)) $(LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_STATIC_LOAD=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER)) $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_SHARED_LOAD=" $(realpath $(PLUMED_SHARED_OBJ)) $(LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_RUNTIME_DEPENDENCIES="   >> $@
	@echo "PLUMED_STATIC_DEPENDENCIES=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))  >> $@
	@echo "PLUMED_SHARED_DEPENDENCIES=" $(realpath $(PLUMED_SHARED_OBJ))  >> $@

Plumed.inc.runtime: Plumed.inc
	@echo Building Plumed.inc, runtime-linking version
	@echo "# PLUMED: runtime installation" > $@
	@cat $< >> $@
	@echo PLUMED_LOAD= '$$(PLUMED_RUNTIME_LOAD)' >> $@
	@echo PLUMED_DEPENDENCIES= '$$(PLUMED_RUNTIME_DEPENDENCIES)' >> $@

Plumed.inc.static: Plumed.inc
	@echo Building Plumed.inc, static version
	@echo "# PLUMED: static installation" > $@
	@cat $< >> $@
	@echo PLUMED_LOAD= '$$(PLUMED_STATIC_LOAD)' >> $@
	@echo PLUMED_DEPENDENCIES= '$$(PLUMED_STATIC_DEPENDENCIES)' >> $@

Plumed.inc.shared: Plumed.inc
	@echo Building Plumed.inc, shared version
	@echo "# PLUMED: shared installation" > $@
	@cat $< >> $@
	@echo PLUMED_LOAD= '$$(PLUMED_SHARED_LOAD)' >> $@
	@echo PLUMED_DEPENDENCIES= '$$(PLUMED_SHARED_DEPENDENCIES)' >> $@

Plumed.cmake: $(OBJ_KERNEL) $(OBJ_WRAPPER) $(PLUMED_SHARED_OBJ) $(OBJ_DYNAMIC_WRAPPER)
	@echo Building Plumed.cmake
	@echo "set(PLUMED_INCLUDE "-I$(realpath ../../include)")" > $@
	@echo "set(PLUMED_RUNTIME_LOAD " $(realpath $(OBJ_DYNAMIC_WRAPPER)) $(LIBS) $(LDFLAGS)")" >> $@
	@echo "set(PLUMED_STATIC_LOAD " $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER)) $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS)")" >> $@
	@echo "set(PLUMED_SHARED_LOAD " $(realpath $(PLUMED_SHARED_OBJ)) $(LIBS) $(LDFLAGS)")" >> $@
	@echo "set(PLUMED_RUNTIME_DEPENDENCIES) "   >> $@
	@echo "set(PLUMED_STATIC_DEPENDENCIES " $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))")"  >> $@
	@echo "set(PLUMED_SHARED_DEPENDENCIES " $(realpath $(PLUMED_SHARED_OBJ))")"  >> $@

Plumed.cmake.runtime: Plumed.cmake
	@echo Building Plumed.cmake, runtime-linking version
	@echo "# PLUMED: runtime installation" > $@
	@cat $< >> $@
	@echo "set("PLUMED_LOAD '$${PLUMED_RUNTIME_LOAD}'")" >> $@
	@echo "set("PLUMED_DEPENDENCIES '$${PLUMED_RUNTIME_DEPENDENCIES}'")" >> $@

Plumed.cmake.static: Plumed.cmake
	@echo Building Plumed.cmake, static version
	@echo "# PLUMED: static installation" > $@
	@cat $< >> $@
	@echo "set("PLUMED_LOAD '$${PLUMED_STATIC_LOAD}'")" >> $@
	@echo "set("PLUMED_DEPENDENCIES '$${PLUMED_STATIC_DEPENDENCIES}'")" >> $@

Plumed.cmake.shared: Plumed.cmake
	@echo Building Plumed.cmake, shared version
	@echo "# PLUMED: shared installation" > $@
	@cat $< >> $@
	@echo "set("PLUMED_LOAD '$${PLUMED_SHARED_LOAD}'")" >> $@
	@echo "set("PLUMED_DEPENDENCIES '$${PLUMED_SHARED_DEPENDENCIES}'")" >> $@


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
	rm -f Plumed.inc* plumed-*

# generic makefile rules
include ../maketools/make.rules


