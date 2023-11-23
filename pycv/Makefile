include Make.inc

#Dependency tracking based on https://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#tldr
#this assumes gcc
DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d

OBJS = PythonCVInterface.o ActionWithPython.o PlumedPythonEmbeddedModule.o

ADDCPPFLAGS=$(shell python3-config --cflags --embed) $(shell python -m pybind11 --includes) $(PLUMED_INCLUDE)
ADDCLDFLAGS=$(shell python3-config --ldflags --embed) 

ifeq ($(SOEXT),dylib)
  SONAME_OPTION:=-Wl,-install_name
else
  SONAME_OPTION:=-Wl,-soname
endif

all: PythonCVInterface.$(SOEXT)


%.o: %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	@echo Compiling object $@
	@$(CXX) -c $(DEPFLAGS) $(CPPFLAGS) $(ADDCPPFLAGS) $(CXXFLAGS) $< -o $@

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(OBJS:%.o=$(DEPDIR)/%.d)
$(DEPFILES):
include $(wildcard $(DEPFILES))

PythonCVInterface.$(SOEXT): $(OBJS)
	@echo Compiling linking $@
	@$(LDSHARED) $(ADDCLDFLAGS) $(SONAME_OPTION),"$(notdir $@)" $(DYNAMIC_LIBS) $(PLUMED_KERNEL) -o $@ $^

clean:
	rm -f $(OBJS) PythonCVInterface.$(SOEXT)