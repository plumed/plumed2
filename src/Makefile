
# include the machine dependent configuration
-include ../Makefile.conf

# if machine dependent configuration has been found:
ifdef GCCDEP

# source for plumed wrapper
SRC_WRAPPER=PlumedStatic.cpp

# source for plumed dynamic wrapper
SRC_DYNAMIC_WRAPPER=Plumed.c

# source for plumed main program
# this should not be linked into libraries
SRC_MAIN=main.cpp

# source for plumed kernel
# everything else
SRC_KERNEL=$(filter-out $(SRC_MAIN) $(SRC_WRAPPER),$(wildcard *.cpp))

# all sources
ALL_SRC=$(SRC_KERNEL) $(SRC_WRAPPER) $(SRC_DYNAMIC_WRAPPER) $(SRC_MAIN)

# corresponding objects:
OBJ_WRAPPER=PlumedStatic.o
OBJ_DYNAMIC_WRAPPER=Plumed.o
OBJ_MAIN=$(SRC_MAIN:.cpp=.o)
OBJ_KERNEL=$(SRC_KERNEL:.cpp=.o)
ALL_OBJ=$(OBJ_KERNEL) $(OBJ_WRAPPER) $(OBJ_DYNAMIC_WRAPPER) $(OBJ_MAIN)

# dependency files for all c++ sources
CPP_SRC=$(wildcard *.cpp)
ALL_DEP=$(addprefix deps/,$(CPP_SRC:.cpp=.d))

# plumed kernel
PLUMED_KERNEL=libplumedKernel.$(SOEXT)

# plumed library (kernel + wrapper)
PLUMED_SHARED_OBJ=libplumed.$(SOEXT)

# plumed main exec
PLUMED_MAIN_STATIC=plumed-static
PLUMED_MAIN_SHARED=plumed
PLUMED_MAIN_RUNTIME=plumed-runtime

# default target
.PHONY: all
all: $(PLUMED_KERNEL) $(PLUMED_SHARED_OBJ) $(OBJ_WRAPPER) \
     $(PLUMED_MAIN_STATIC) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME) \
     Plumed.inc Plumed.inc.runtime Plumed.inc.shared Plumed.inc.static \
     compile_options.sh

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
	@echo "PLUMED_INCLUDE="-I$(realpath ../include) > $@
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

$(ALL_DEP): PlumedConfig.h Makefile.conf.xxd

# this compiler is used for dependencies generation
# if unset, use normal CXX
ifndef CXXDEP
CXXDEP=$(CXX)
endif

# rule for c++
%.o: %.cpp
	$(CXXDEP) -c -MM -MF$*.d $(CPPFLAGS) $(CXXFLAGS) $*.cpp
	@cp -f $*.d $*.d.tmp
	@-test -d deps || mkdir deps
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	 sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@mv $*.d deps/$*.d
	@rm -f $*.d.tmp
ifndef XLF
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $*.cpp -o $*.o
endif

# rule for c
%.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $*.c -o $*.o

# there is a single C file here, and this is its dependency:
Plumed.o: Plumed.c Plumed.h

# file to import compilation options inside c++
PlumedConfig.h: PlumedConfig.h.in ../Makefile.conf
	@echo Building PlumedConfig.h, an includable c++ header with some compilation options
	@sed "s/@SOEXT@/$(SOEXT)/g" $< | \
         sed "s|@PLUMED_ROOT@|$(PWD:src=)|g" > $@

# file to import compilation options inside a bash script
compile_options.sh:
	@echo Building compile_options.sh, a sourceable bash script with some compilation options
	@echo "compile='"$(CXX) -c -I$(realpath .) $(CPPFLAGS) $(CXXFLAGS) -o"'" > $@
	@echo "link='"$(LDSO) -o"'" >> $@
	@echo "soext='"$(SOEXT)"'" >> $@

Makefile.conf.xxd: ../Makefile.conf
	cat ../Makefile.conf | xxd -i > Makefile.conf.xxd

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
	rm -f *.o
	rm -f *.d
	rm -fr deps
	rm -f *.cdep
	rm -f *.cppdep
	rm -f PlumedConfig.h
	rm -f Plumed.inc*
	rm -f *~
	rm -f *.tmp
	rm -f *.so
	rm -f *.dylib
	rm -f *.xxd
	rm -f compile_options.sh
	rm -f $(PLUMED_MAIN_STATIC) $(PLUMED_MAIN_SHARED) $(PLUMED_MAIN_RUNTIME)

