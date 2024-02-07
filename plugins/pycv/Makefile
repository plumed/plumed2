# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../../Makefile.conf
	ifndef canPyCV
  	include ./Makefile.conf
	endif
endif

.PHONY: clean check all
#Dependency tracking based on https://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr
#this assumes gcc
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
PYBIND11FLAG = 
ADDCPPFLAGS:=$(python_cf_embedded) $(pybind11_cflags) $(ADDCPPFLAGS)
ADDCLDFLAGS:=$(ADDCLDFLAGS) $(python_ld_embedded)
OBJS = ActionWithPython.o PythonCVInterface.o PythonFunction.o PlumedPythonEmbeddedModule.o

ifeq ($(SOEXT),dylib)
  SONAME_OPTION:=-Wl,-install_name
else
  SONAME_OPTION:=-Wl,-soname
endif

all: PythonCVInterface.$(SOEXT)

#-fvisibility=hidden is needed for pybind11 (to not conflict with different pybind11 versions)
#I think I enforced this nearly everywhere I set up a flag for the compiler
ActionWithPython.o PythonCVInterface.o PythonFunction.o PlumedPythonEmbeddedModule.o: PYBIND11FLAG= -fvisibility=hidden

%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling object $@
	@$(CXX) -c $(DEPFLAGS) $(CPPFLAGS) $(PYBIND11FLAG) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@


$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(OBJS:%.o=$(DEPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

#-Wl,--no-as-needed forces the python library to be linked, without this in a WSL does not work
#TODO: seems that $PLUMED_KERNEL is not needed, check
PythonCVInterface.$(SOEXT): $(OBJS)
	@echo Linking $@
	$(LDSHARED) $(SONAME_OPTION),"$(notdir $@)" $(DYNAMIC_LIBS) $(PLUMED_KERNEL) $(ADDCLDFLAGS) $^ -o $@ 

clean:
	rm -f $(OBJS) PythonCVInterface.$(SOEXT)

check: all
	$(MAKE) -C regtest testclean
	$(MAKE) -C regtest checkfail
