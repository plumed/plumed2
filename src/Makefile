
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

# dependency files
ALL_DEP=$(ALL_OBJ:.o=.d)

# plumed kernel
PLUMED_KERNEL=libplumedKernel.$(SOEXT)

# plumed library (kernel + wrapper)
PLUMED_SHARED_OBJ=libplumed.$(SOEXT)

# plumed main exec
PLUMED_MAIN_STATIC=plumed-static
PLUMED_MAIN_SHARED=plumed

# default target
.PHONY: all
all: $(PLUMED_KERNEL) $(PLUMED_SHARED_OBJ) $(OBJ_WRAPPER) \
     $(PLUMED_MAIN_STATIC) $(PLUMED_MAIN_SHARED) \
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

Plumed.inc: $(OBJ_KERNEL) $(OBJ_WRAPPER) $(PLUMED_SHARED_OBJ) $(OBJ_DYNAMIC_WRAPPER)
	@echo Building Plumed.inc
	@echo "PLUMED_INCLUDE="-I$(realpath ../include) >> $@
	@echo "PLUMED_RUNTIME_LOAD=" $(realpath $(OBJ_DYNAMIC_WRAPPER)) $(LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_STATIC_LOAD=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER)) $(LIBS) $(DYNAMIC_LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_SHARED_LOAD=" $(realpath $(PLUMED_SHARED_OBJ)) $(LIBS) $(LDFLAGS) >> $@
	@echo "PLUMED_RUNTIME_DEPENDENCIES="   >> $@
	@echo "PLUMED_STATIC_DEPENDENCIES=" $(realpath $(OBJ_KERNEL)) $(realpath $(OBJ_WRAPPER))  >> $@
	@echo "PLUMED_RUNTIME_DEPENDENCIES=" $(realpath $(PLUMED_SHARED_OBJ))  >> $@

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

$(ALL_DEP): PlumedConfig.h

# rule for c++
%.o: %.cpp
	$(CXX) -c -MM -MF$*.d $(CPPFLAGS) $(CXXFLAGS) $*.cpp
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	 sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp
ifndef XLF
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $*.cpp -o $*.o
endif

# rule for c
%.o: %.c
	$(CC) -c -MM -MF$*.d -M $(CPPFLAGS) $(CFLAGS) $*.c
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	 sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp
ifndef XLF
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $*.c -o $*.o
endif

# file to import compilation options inside c++
PlumedConfig.h: PlumedConfig.h.in ../Makefile.conf
	@echo Building PlumedConfig.h, an includable c++ header with some compilation options
	@sed "s/@SOEXT@/$(SOEXT)/g" $< > $@

# file to import compilation options inside a bash script
compile_options.sh:
	@echo Building compile_options.sh, a sourceable bash script with some compilation options
	@echo "compile='"$(CXX) -c -I$(realpath .) $(CPPFLAGS) $(CXXFLAGS) -o"'" > $@
	@echo "link='"$(LDSO) -o"'" >> $@
	@echo "soext='"$(SOEXT)"'" >> $@

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
	rm -f *.cdep
	rm -f *.cppdep
	rm -f PlumedConfig.h
	rm -f Plumed.inc*
	rm -f *~
	rm -f *.tmp
	rm -f *.so
	rm -f *.dylib
	rm -f compile_options.sh
	rm -f $(PLUMED_MAIN_STATIC) $(PLUMED_MAIN_SHARED)

