include Make.inc

#Dependency tracking based on https://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr
#this assumes gcc
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
PYBIND11FLAG = 

OBJS = ActionWithPython.o PythonCVInterface.o PythonFunction.o PlumedPythonEmbeddedModule.o

ADDCPPFLAGS=$(shell python3-config --cflags --embed) $(shell python3 -m pybind11 --includes) $(PLUMED_INCLUDE)
ADDCLDFLAGS=$(shell python3-config --ldflags --embed) 

ifeq ($(SOEXT),dylib)
  SONAME_OPTION:=-Wl,-install_name
else
  SONAME_OPTION:=-Wl,-soname
endif

all: PythonCVInterface.$(SOEXT)

#-fvisibility=hidden is needed for pybind11 (to not conflict with different pybind11 versions)
ActionWithPython.o PythonCVInterface.o PythonFunction.o PlumedPythonEmbeddedModule.o: PYBIND11FLAG= -fvisibility=hidden

%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling object $@
	@$(CXX) -c $(DEPFLAGS) $(CPPFLAGS) $(PYBIND11FLAG) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@


$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(OBJS:%.o=$(DEPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

PythonCVInterface.$(SOEXT): $(OBJS)
	@echo Compiling linking $@
	@$(LDSHARED) $(ADDCLDFLAGS) $(SONAME_OPTION),"$(notdir $@)" $(DYNAMIC_LIBS) $(PLUMED_KERNEL) -o $@ $^

clean:
	rm -f $(OBJS) PythonCVInterface.$(SOEXT)

check: all
	$(MAKE) -C regtest
